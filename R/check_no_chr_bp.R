#' Ensure that CHR and BP are missing if SNP is present, can find them
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return A list containing two data tables:
#' \describe{
#'   \item{\code{sumstats_dt}}{The modified summary statistics data table.}
#'   \item{\code{rsids}}{`snpsById`, filtered to SNPs of interest if loaded 
#'   already, or else `NULL`.}
#'   \item{\code{log_files}}{List of log files.}
#' }
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table setcolorder
#' @importFrom BSgenome snpsById
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_no_chr_bp <- function(sumstats_dt, 
                            path, 
                            ref_genome,
                            rsids,
                            imputation_ind,
                            log_folder_ind,
                            check_save_out,
                            tabix_index, 
                            nThread,
                            log_files,
                            dbSNP,
                            dbSNP_tarball) {
    SNP <- i.seqnames <- CHR <- BP <- i.pos <- LP <- P <- IMPUTATION_CHR <-
        IMPUTATION_BP <- NULL
    # If SNP present but no CHR/BP then need to find them
    col_headers <- names(sumstats_dt)
    if (sum(c("CHR", "BP") %in% col_headers) <= 1 &
        sum("SNP" %in% col_headers) == 1) {
        # if dataset has one of CHR or BP remove it and take from re dataset
        if (sum(c("CHR", "BP") %in% col_headers) == 1) {
            colsToDelete <- c("CHR", "BP")[c("CHR", "BP") %in% col_headers]
            sumstats_dt[, (colsToDelete) := NULL]
        }
        # check if rsids loaded if not do so
        if (is.null(rsids)) {
            rsids <-
                load_ref_genome_data(
                    data.table::copy(sumstats_dt$SNP), 
                    ref_genome = ref_genome,
                    dbSNP = dbSNP,
                    dbSNP_tarball = dbSNP_tarball,
                    msg="Chromosome or Base Pair Position"
                )
        } else {
            print_msg <- paste0(
                "There is no Chromosome or Base Pair Position column",
                " found within the data. It must be inferred from ",
                " other column information."
            )
            message(print_msg)
        }
        # ensure rsids is up-to-date with filtered sumstats_dt
        rsids <- rsids[unique(sumstats_dt$SNP), , nomatch = NULL]
        data.table::setkey(rsids, SNP)
        # join on CHR BP to sumstats
        data.table::setkey(sumstats_dt, SNP)
        sumstats_dt[rsids, CHR := i.seqnames]
        sumstats_dt[rsids, BP := i.pos]
        # remove rows where CHR/BP couldn't be found
        # If user wants log, save it to there
        if (log_folder_ind) {
            name <- "chr_bp_not_found_from_snp"
            name <- get_unique_name_log_file(name = name,
                                             log_files = log_files)
            write_sumstats(
                sumstats_dt =
                    sumstats_dt[!complete.cases(sumstats_dt[, c(
                        "CHR",
                        "BP"
                    )]), ],
                save_path =
                    paste0(
                        check_save_out$log_folder,
                        "/", name,
                        check_save_out$extension
                    ),
                sep = check_save_out$sep,
                #don't tab indx as could be miss values & cause err
                #tabix_index = tabix_index,
                nThread = nThread
            )
            log_files[[name]] <-
                paste0(check_save_out$log_folder, "/",
                       name, check_save_out$extension)
        }
        sumstats_dt <- sumstats_dt[complete.cases(
            sumstats_dt[, c("CHR", "BP")]), ]
        # move SNP, CHR, BP to start
        other_cols <-
            names(sumstats_dt)[!names(sumstats_dt) %in% c("SNP", "CHR", "BP")]
        data.table::setcolorder(sumstats_dt,
                                c("SNP", "CHR", "BP", other_cols))

        # if user specifies add a column to notify of the imputation
        if (imputation_ind) {
            sumstats_dt[, IMPUTATION_CHR := TRUE]
            sumstats_dt[, IMPUTATION_BP := TRUE]
        }

        return(list("sumstats_dt" = sumstats_dt,
                    "rsids" = rsids, 
                    "log_files" = log_files))
    }
    #also check if bp and chr are present but values are na
    if(sum("SNP" %in% col_headers) == 1){
      na_chr_bp <- 
        sumstats_dt[grep("^rs", SNP),][is.na(CHR) & is.na(BP),]
      if (nrow(na_chr_bp)>0 && 
          sum(c("CHR", "BP","SNP") %in% col_headers) == 3){
        #impute the chr,bp data
        # if dataset has one of CHR or BP remove it and take from re dataset
        if (sum(c("CHR", "BP") %in% col_headers) >= 1) {
          colsToDelete <- c("CHR", "BP")[c("CHR", "BP") %in% col_headers]
          na_chr_bp[, (colsToDelete) := NULL]
        }
        #pass the dataset without chr bp to be imputed
        #keep old rsids to update
        old_rsids <- na_chr_bp$SNP
        na_chr_bp <-
          check_no_chr_bp(sumstats_dt=na_chr_bp,path=path,ref_genome=ref_genome,
                          rsids=rsids,imputation_ind=imputation_ind,
                          log_folder_ind=log_folder_ind,
                          check_save_out=check_save_out,
                          tabix_index=tabix_index,nThread=nThread,
                          log_files=log_files,dbSNP=dbSNP,
                          dbSNP_tarball=dbSNP_tarball) 
        #join back on the sumstats
        sumstats_dt <- sumstats_dt[!SNP %in% na_chr_bp$sumstats_dt$SNP,]
        sumstats_dt <- rbindlist(list(sumstats_dt,na_chr_bp$sumstats_dt),
                                 fill=TRUE)
        #join on the log files
        log_files <- c(log_files,na_chr_bp$log_files)
        #update rsids
        old_rsids_drop <- old_rsids[!old_rsids %in% na_chr_bp$sumstats_dt]
        rsids <- rsids[!rsids %in% old_rsids_drop]
      }
    }
    return(list("sumstats_dt" = sumstats_dt, 
                "rsids" = rsids,
                "log_files" = log_files))
}
