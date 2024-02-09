#' Ensure that SNP appears to be valid RSIDs (starts with rs)
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list.
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table :=
#' @importFrom data.table setcolorder
#' @importFrom data.table copy
#' @importFrom data.table tstrsplit
#' @importFrom data.table rbindlist
#' @importFrom BSgenome snpsByOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom stringr str_sub
check_no_rs_snp <- function(sumstats_dt, path, ref_genome, snp_ids_are_rs_ids,
                            indels,imputation_ind, log_folder_ind, 
                            check_save_out,tabix_index, nThread, log_files,
                            dbSNP) {
    SNP <- CHR <- CHR1 <- BP1 <- i.RefSNP_id <- IMPUTATION_SNP <-
        SNP_old_temp <- SNP_INFO <- A1 <- A2 <- .I <- NULL
    # if snp ids aren't rs ids rename the column to ID's 
    # so RSIDs can be inferred
    if ((!snp_ids_are_rs_ids) & sum("SNP" %in% names(sumstats_dt)) == 1) {
        data.table::setnames(sumstats_dt, "SNP", "ID")
    }
    #also check if user didn't set the param snp_ids_are_rs_ids but it's clear
    if("SNP" %in% names(sumstats_dt)){
      #If SNP present, make sure it isn't a proxy for CHR:bp
      #RS ID should be all numeric or should start with rs and then be numeric
      snps <- sumstats_dt$SNP
      snp_check <- suppressMessages(as.numeric(
        stringr::str_sub(string = snps, start = 3, end = -1L) 
      ))
      #this won't catch chr:bp when chr is just number i.e. 1:123456789
      #so second check needed
      snp_check2 <- unique(c(grep(":", snps),
                             grep("-", snps),
                             grep("_", snps)))
      #however if they start with rs they stay i.e. rs140052487:C:A
      snp_check3 <- grep("^rs", snps)
      #remove these from those that passed first check
      if(length(snp_check2)!=0)
        snp_check <- snp_check[-snp_check2]
      num_rsids <- length(snps[!is.na(snp_check)])
      #if num_rsids==0 then they clearly meant SNP col to not represent 
      #RSIDs so rename to avoid errors later
      if(num_rsids==0 && length(snp_check3)==0)
        setnames(sumstats_dt,"SNP","ID")
    }
    # If SNP column doesn't start with rs
    col_headers <- names(sumstats_dt)
    if (sum("SNP" %in% col_headers) == 1) {
        message("Checking SNP RSIDs.")
        # needed for later to join and match SNPs
        if (imputation_ind) {
            # if any are NA, we need to give them a unique ID
            # need to be able to identify to revert back so picked up 
            #by other, NA for RS ID check
            sumstats_dt[is.na(SNP)|SNP=="", SNP := paste0("NA_",.I)]
            sumstats_dt[, SNP_old_temp := SNP]
        }
        miss_rs <- sumstats_dt[!grep("^rs", SNP), ]
        # first case is something other than rs id and chr:bp - impute SNP
        # second case is chr:bp together - impute SNP for these
        miss_rs_chr_bp <- miss_rs[grep(":", SNP),] 
        miss_rs_chr_bp <- miss_rs_chr_bp[!grep(".*:.*:.*", SNP),]
        if (nrow(miss_rs) != nrow(sumstats_dt) && nrow(miss_rs_chr_bp) > 0) {
            # check if impute of correct SNP ID possible
            if (sum(c("CHR", "BP") %in% col_headers) == 2 &&
                nrow(miss_rs[!grep(":", SNP), ]) > 0) {
                bad_snp <- miss_rs[!grep(":", SNP), ]
                msg <- paste0(
                    formatC(nrow(bad_snp), big.mark = ","),
                    " SNP IDs are not correctly formatted.",
                    " These will be corrected from the reference genome."
                )
                message(msg)
                # remove snp column and pass to function to impute snp
                bad_snp <- bad_snp[, SNP := NULL]
                # now impute correct RS ID for those missing it
                corrected_snp <-
                    check_no_snp(
                        sumstats_dt = bad_snp, path = tempfile(),
                        ref_genome = ref_genome,
                        indels = indels,
                        imputation_ind = imputation_ind,
                        log_folder_ind = log_folder_ind,
                        check_save_out = check_save_out, 
                        tabix_index = tabix_index,
                        nThread = nThread, 
                        log_files = log_files, 
                        dbSNP = dbSNP,
                        verbose = FALSE
                    )
                log_files <- corrected_snp$log_files
                corrected_snp <- corrected_snp$sumstats_dt
                # make sure columns in correct order
                data.table::setcolorder(corrected_snp, names(sumstats_dt))
                # remove rows missing from the reference genome and combine
                # If IMPUTATION column added add it to other DT
                if (imputation_ind &&
                    !"IMPUTATION_SNP" %in% names(sumstats_dt)) {
                    sumstats_dt[, IMPUTATION_SNP := NA]
                }
                sumstats_dt <-
                    data.table::rbindlist(list(
                        sumstats_dt[grep("^rs", SNP), ],
                        corrected_snp
                    ))
            } else {
                # remove snps missing rs
                # If user wants log, save it to there
                if (log_folder_ind) {
                    name <- "snp_missing_rs"
                    name <- get_unique_name_log_file(name = name,
                                                     log_files = log_files)
                    write_sumstats(
                        sumstats_dt = sumstats_dt[!grepl("^rs", SNP), ],
                        save_path =
                            paste0(
                                check_save_out$log_folder,
                                "/", name,
                                check_save_out$extension
                            ),
                        sep = check_save_out$sep,
                        tabix_index = tabix_index,
                        nThread = nThread
                    )
                    log_files[[name]] <-
                        paste0(check_save_out$log_folder,
                               "/", name, check_save_out$extension)
                }
                sumstats_dt <- sumstats_dt[grep("^rs", SNP), ]
            }
            # check if any have more than 1 ":" remove these
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "snp_multi_colon"
                name <- get_unique_name_log_file(name = name,
                                                 log_files = log_files)
                write_sumstats(
                    sumstats_dt = miss_rs_chr_bp[grep(".*:.*:.*", SNP)],
                    save_path =
                        paste0(
                            check_save_out$log_folder,
                            "/", name,
                            check_save_out$extension
                        ),
                    sep = check_save_out$sep,
                    tabix_index = tabix_index,
                    nThread = nThread
                )
                log_files[[name]] <-
                    paste0(
                        check_save_out$log_folder, "/", name,
                        check_save_out$extension
                    )
            }
            #remove cases with more than one colon i.e. not just chr:bp
            miss_rs_chr_bp <- miss_rs_chr_bp[!grep(".*:.*:.*", SNP)]
            msg <- paste0(
                formatC(nrow(miss_rs_chr_bp), big.mark = ","),
                " SNP IDs appear to be made up of ",
                "chr:bp, these will be replaced by their SNP ID from the",
                " reference genome"
            )
            message(msg)
            SNP_LOC_DATA <- load_snp_loc_data(ref_genome,dbSNP, NULL)
            # split out chr:bp - check if chr or bp first by longer of two
            splits <- strsplit(miss_rs_chr_bp[1, SNP],
                               split = ":", fixed = TRUE)[[1]]
            if (nchar(splits[1]) < nchar(splits[2])) {
                format <- c("CHR1", "BP1")
            } else {
                format <- c("BP1", "CHR1")
            }
            miss_rs_chr_bp[, (format) := data.table::tstrsplit(SNP,
                                                               split = ":", 
                                                               fixed = TRUE)]
            # if BP col has other info after, drop it
            if (sum(grepl("[[:punct:]].*", miss_rs_chr_bp$BP1)) > 0) {
                miss_rs_chr_bp[, BP1 := gsub("([[:punct:]]).*", "", BP1)]
            }
            # ensure integer col
            miss_rs_chr_bp[, BP1 := as.integer(BP1)]
            # now drop SNP
            miss_rs_chr_bp[, SNP := NULL]
            # if chromosome col has chr prefix remove it
            miss_rs_chr_bp[, CHR1 := gsub("chr", "", CHR1)]
            # if chromosome col has other info after, drop it
            if (sum(grepl("[[:punct:]].*", miss_rs_chr_bp$CHR1)) > 0) {
                miss_rs_chr_bp[, CHR1 := gsub("([[:punct:]]).*", "", CHR1)]
            }
            # avoid SNPs with NA values in chr or bp
            gr_snp <- data.table::copy(miss_rs_chr_bp)
            # if indels are in datqset, these should not be imputed as RS ID
            # will instead relate to a SNP at same position - remove these
            if(sum(c("A1", "A2") %in% col_headers) == 2 & indels ){
                #identify Indels based on num char in A1, A2
                num_indels <- nrow(gr_snp[(nchar(A1)>1 | nchar(A2)>1),])
                if(num_indels>0){
                    msg <- paste0("Found ",
                                  formatC(nrow(num_indels),big.mark = ","),
                                  " Indels. These won't",
                                  " be checked against the reference ",
                                  "genome as it does not contain ",
                                  "Indels.\nWARNING If your sumstat ",
                                  "doesn't contain Indels, set the ",
                                  "indel param to FALSE & rerun ",
                                  "MungeSumstats::format_sumstats()")
                    message(msg)
                    gr_snp <- gr_snp[!(nchar(A1)>1 | nchar(A2)>1),]
                }
            }
            if(nrow(gr_snp)>0){
                incl_cols <- c("CHR1", "BP1")
                gr_snp <- gr_snp[complete.cases(
                    gr_snp[, incl_cols, with = FALSE])]
                gr_snp <-
                    GenomicRanges::makeGRangesFromDataFrame(gr_snp,
                                                            keep.extra.columns = TRUE,
                                                            seqnames.field = "CHR1",
                                                            start.field = "BP1",
                                                            end.field = "BP1"
                    )
                gr_rsids <-
                    BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = gr_snp)
                rsids <- data.table::setDT(data.frame(gr_rsids))
                data.table::setnames(rsids, "seqnames", "CHR1")
                data.table::setnames(rsids, "pos", "BP1")
                # in case there is CHR8 and chr8
                rsids[, CHR1 := tolower(as.character(CHR1))]
                miss_rs_chr_bp[, CHR1 := tolower(as.character(CHR1))]
                # join on SNP ID to sumstats
                data.table::setkeyv(miss_rs_chr_bp, c("CHR1", "BP1"))
                data.table::setkeyv(rsids, c("CHR1", "BP1"))
                miss_rs_chr_bp[rsids, SNP := i.RefSNP_id]
                # remove rows where SNP couldn't be found
                # If user wants log, save it to there
                if (log_folder_ind) {
                    name <- "snp_not_found_from_bp_chr"
                    name <- get_unique_name_log_file(name = name,
                                                     log_files = log_files)
                    write_sumstats(
                        sumstats_dt =
                            miss_rs_chr_bp[!complete.cases(
                                miss_rs_chr_bp[, "SNP"]), ],
                        save_path =
                            paste0(
                                check_save_out$log_folder,
                                "/", name,
                                check_save_out$extension
                            ),
                        sep = check_save_out$sep,
                        tabix_index = tabix_index,
                        nThread = nThread
                    )
                    log_files[[name]] <-
                        paste0(
                            check_save_out$log_folder, "/", name,
                            check_save_out$extension
                        )
                }
                miss_rs_chr_bp <- miss_rs_chr_bp[complete.cases(
                    miss_rs_chr_bp[, "SNP"]), ]
                # remove temp columns
                miss_rs_chr_bp[, (format) := NULL]
                # join with full dataset
                # If IMPUTATION column added add it to other DT
                if (imputation_ind && 
                    "IMPUTATION_SNP" %in% names(sumstats_dt) &&
                    !"IMPUTATION_SNP" %in% names(miss_rs_chr_bp)) {
                    miss_rs_chr_bp[, IMPUTATION_SNP := NA]
                }
                # get columns in same order as rest of data table
                data.table::setcolorder(miss_rs_chr_bp, col_headers)
                sumstats_dt <- data.table::rbindlist(
                    list(sumstats_dt, miss_rs_chr_bp))
            }    
        }
        if (nrow(miss_rs) != nrow(sumstats_dt) && nrow(miss_rs) != 0) {
            #Prev if (nrow(miss_rs_chr_bp) == 0) {
            if (nrow(miss_rs_chr_bp) < nrow(miss_rs)) {
                # don't filter twice if hit prev condition
                # check if impute of correct SNP ID possible
                if (sum(c("CHR", "BP") %in% col_headers) == 2 &&
                    nrow(sumstats_dt[!grep("^rs", SNP), ]) > 0) {
                    bad_snp <- sumstats_dt[!grep("^rs", SNP), ]
                    msg <- paste0(
                        formatC(nrow(bad_snp), big.mark = ","),
                        " SNP IDs are not correctly formatted.",
                        " These will be corrected from the reference genome."
                    )
                    message(msg)
                    # remove snp column and pass to function to impute snp
                    bad_snp <- bad_snp[, SNP := NULL]
                    if(sum(c("A1", "A2") %in% col_headers) == 2 & indels ){
                        #identify Indels based on num char in A1, A2
                        num_indels <- nrow(bad_snp[(nchar(A1)>1 | 
                                                        nchar(A2)>1),])
                        if(num_indels>0){
                            msg <- paste0("Found ",
                                          formatC(nrow(num_indels),big.mark = ","),
                                          " Indels. These won't",
                                          " be checked against the referenc",
                                          "e genome as it does not contain ",
                                          "Indels.\nWARNING If your sumstat ",
                                          "doesn't contain Indels, set the ",
                                          "indel param to FALSE & rerun ",
                                          "MungeSumstats::format_sumstats()")
                            message(msg)
                            bad_snp <- bad_snp[!(nchar(A1)>1 | nchar(A2)>1),]
                        }
                    }
                    # now impute correct RS ID for those missing it
                    corrected_snp <-
                        check_no_snp(
                            sumstats_dt = bad_snp, path = tempfile(),
                            ref_genome = ref_genome, 
                            indels = indels,
                            imputation_ind = imputation_ind,
                            log_folder_ind = log_folder_ind,
                            check_save_out = check_save_out,
                            tabix_index = tabix_index,
                            nThread = nThread, 
                            log_files = log_files,
                            dbSNP=dbSNP,
                            verbose = FALSE
                        )
                    log_files <- corrected_snp$log_files
                    corrected_snp <- corrected_snp$sumstats_dt
                    # make sure columns in correct order
                    data.table::setcolorder(corrected_snp,
                                            names(sumstats_dt))
                    # remove rows missing from the reference genome and combine
                    # If IMPUTATION column added add it to other DT
                    if (imputation_ind &&
                        !"IMPUTATION_SNP" %in% names(sumstats_dt)) {
                        sumstats_dt[, IMPUTATION_SNP := NA]
                    }
                    sumstats_dt <-
                        data.table::rbindlist(list(
                            sumstats_dt[grep("^rs", SNP), ],
                            corrected_snp
                        ))
                }
                # remove snps missing rs
                # If user wants log, save it to there
                if (log_folder_ind &&
                    nrow(sumstats_dt[!grep("^rs", SNP), ]) > 0) {
                    #if imputation ind on, NA's will be replaced with NA_i
                    #change back here
                    tmp_rmv <- sumstats_dt[!grepl("^rs", SNP), ]
                    if(imputation_ind){
                      tmp_rmv[, SNP_old_temp := NULL]
                      tmp_rmv[grepl("^NA_",SNP), SNP := NA]
                    }
                    name <- "snp_missing_rs"
                    name <- get_unique_name_log_file(name = name,
                                                     log_files = log_files)
                    write_sumstats(
                        sumstats_dt = tmp_rmv,
                        save_path =
                            paste0(
                                check_save_out$log_folder,
                                "/", name,
                                check_save_out$extension
                            ),
                        sep = check_save_out$sep,
                        tabix_index = tabix_index,
                        nThread = nThread
                    )
                    log_files[[name]] <-
                        paste0(check_save_out$log_folder, "/",
                               name, check_save_out$extension)
                }
                sumstats_dt <- sumstats_dt[grep("^rs", SNP), ]
            }
            # if any weird SNP rows left that aren't 
            # chr:bp or rs id's remove them
            if (sum(c("CHR", "BP") %in% col_headers) != 2) {
                msg <-
                    paste0(
                        formatC(nrow(miss_rs) - nrow(miss_rs_chr_bp),
                                big.mark = ","),
                        " SNP IDs are not correctly",
                        " formatted and will be removed."
                    )
                message(msg)
            }
        }
        # also check for weirdly formatted SNPs like rs1234:.....
        # remove everything after : and put in separate column
        # (can extend code to infer info from this separate column later)
        rs_plus <- sumstats_dt[intersect(grep("^rs", SNP), grep(":", SNP)), ]
        if (nrow(rs_plus) != 0) {
            msg <- paste0(
                formatC(nrow(rs_plus), big.mark = ","),
                " SNP IDs contain other information in the same column.",
                " These will be separated."
            )
            message(msg)
            setnames(rs_plus, "SNP", "SNP_INFO")
            # isolate RS ID
            rs_plus[, SNP := sub("\\:.*", "", SNP_INFO)]
            # isolate extra info
            rs_plus[, SNP_INFO := sub("^.*?:", "", SNP_INFO)]
            # join back
            # add extra column
            sumstats_dt[, SNP_INFO := NA]
            # make sure columns in correct order
            data.table::setcolorder(rs_plus, names(sumstats_dt))
            # update values
            sumstats_dt <-
                data.table::rbindlist(list(
                    sumstats_dt[!intersect(
                        grep("^rs", SNP),
                        grep(":", SNP)
                    ), ],
                    rs_plus
                ))
        }
        # if imputation_ind return column specifying imputed
        if (imputation_ind) {
            # make sure there are rows that would need imputing
            if (nrow(miss_rs) > 0 && nrow(miss_rs) != nrow(sumstats_dt) &&
                (nrow(miss_rs_chr_bp) > 0 | nrow(miss_rs) != 0)) {
                setkey(miss_rs, SNP_old_temp)
                setkey(sumstats_dt, SNP_old_temp)
                sumstats_dt[miss_rs, IMPUTATION_SNP := TRUE]
            }
            # remove temp columns either way
            sumstats_dt[, SNP_old_temp := NULL]
        }
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}