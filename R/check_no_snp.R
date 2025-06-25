#' Ensure that SNP is present if not can find it with CHR and BP
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @param verbose should messages be printed. Default it TRUE.
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log files list
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table :=
#' @importFrom data.table setcolorder
#' @importFrom data.table copy
#' @importFrom BSgenome snpsByOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
check_no_snp <- function(sumstats_dt, path, ref_genome, snp_ids_are_rs_ids, 
                         indels, imputation_ind, log_folder_ind, check_save_out,
                         tabix_index, nThread, log_files, dbSNP,
                         dbSNP_tarball = NULL,msg = NULL,verbose = TRUE) {
    SNP <- CHR <- i.RefSNP_id <- IMPUTATION_SNP <- BP <- A1 <- A2 <- NULL
    # If CHR and BP are present BUT not SNP then need 
    # to find the relevant SNP ids
    col_headers <- names(sumstats_dt)
    if (sum(c("CHR", "BP") %in% col_headers) == 2 &
        sum("SNP" %in% col_headers) == 0) {
        msg <- "SNP"
        if (isFALSE(verbose)) {
            msg <- NULL
        }
        SNP_LOC_DATA <- load_snp_loc_data(
          ref_genome = ref_genome,
          dbSNP = dbSNP,
          dbSNP_tarball = dbSNP_tarball,
          msg = NULL
        )
        # if chromosome col has chr prefix remove it
        sumstats_dt[, CHR := gsub("chr", "", CHR)]
        # avoid SNPs with NA values in chr or bp
        gr_snp <- data.table::copy(sumstats_dt)
        # if indels are in datqset, these should not be imputed as RS ID
        # will instead relate to a SNP at same position - remove these
        if(sum(c("A1", "A2") %in% col_headers) == 2 & indels ){
            #identify Indels based on num char in A1, A2
            num_indels <- nrow(gr_snp[(nchar(A1)>1 | nchar(A2)>1),])
            if(num_indels>0){
                msg <- paste0("Found ",
                          formatC(num_indels,big.mark = ","),
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
        incl_cols <- c("CHR", "BP")
        gr_snp <- gr_snp[complete.cases(gr_snp[, incl_cols, with = FALSE])]
        #could be 0 with indels, causes error
        if(nrow(gr_snp)>0){
            gr_snp <-
                GenomicRanges::makeGRangesFromDataFrame(gr_snp,
                    keep.extra.columns = TRUE,
                    seqnames.field = "CHR",
                    start.field = "BP",
                    end.field = "BP"
                )
            gr_rsids <-
                BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = gr_snp)
            rsids <- data.table::setDT(data.frame(gr_rsids))
        }
        else{
            #just add fake rsids so no match found
            rsids <-
                data.table::data.table(seqnames=0,pos=0,strand="*",
                                       RefSNP_id="fake",alleles_as_ambig="Y")
        }
        data.table::setnames(rsids, "seqnames", "CHR")
        data.table::setnames(rsids, "pos", "BP")
        # in case there is CHR8 and chr8 - but keep sex chr as upper
        rsids[, CHR := tolower(as.character(CHR))]
        rsids[, CHR := gsub("x|23", "X", CHR)]
        rsids[, CHR := gsub("y", "Y", CHR)]
        sumstats_dt[, CHR := tolower(as.character(CHR))]
        sumstats_dt[, CHR := gsub("x|23", "X", CHR)]
        sumstats_dt[, CHR := gsub("y", "Y", CHR)]
        # ensure bp is numeric
        sumstats_dt[, BP := as.numeric(BP)]
        # join on SCHR BP to sumstats
        data.table::setkeyv(sumstats_dt, c("CHR", "BP"))
        data.table::setkeyv(rsids, c("CHR", "BP"))
        sumstats_dt[rsids, SNP := i.RefSNP_id]
        # remove rows where SNP couldn't be found
        # If user wants log, save it to there
        if (log_folder_ind) {
            name <- "snp_not_found_from_chr_bp"
            name <- get_unique_name_log_file(name = name,
                                             log_files = log_files)
            write_sumstats(
                sumstats_dt =
                    sumstats_dt[!complete.cases(sumstats_dt[, c("SNP")]), ],
                save_path =
                    paste0(
                        check_save_out$log_folder,
                        "/", name,
                        check_save_out$extension
                    ),
                sep = check_save_out$sep,
                #don't tab index these as could be missing values and cause err
                #tabix_index = tabix_index,
                nThread = nThread
            )
            log_files[[name]] <-
                paste0(check_save_out$log_folder, "/",
                       name, check_save_out$extension)
        }
        sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt[, "SNP"]), ]
        # move SNP to start
        other_cols <- names(sumstats_dt)[names(sumstats_dt) != "SNP"]
        data.table::setcolorder(sumstats_dt, c("SNP", other_cols))

        # if user specifies add a column to notify of the imputation
        if (imputation_ind) {
            sumstats_dt[, IMPUTATION_SNP := TRUE]
        }
        return(list("sumstats_dt" = sumstats_dt,
                    "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt,
                    "log_files" = log_files))
    }
}
