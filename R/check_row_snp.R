#' Ensure all rows have SNPs beginning with rs or SNP, drop those that don't
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and log file list
#' @keywords internal
#' @importFrom data.table fread
#' @importFrom data.table fwrite
check_row_snp <- function(sumstats_dt, path, log_folder_ind, check_save_out,
                          tabix_index, nThread, log_files) {
    SNP <- NULL
    # All rows should start with either SNP or rs... if they don't drop them
    # use data table for speed
    num_bad_ids <- nrow(sumstats_dt[!grep("^rs", SNP), ])
    if (num_bad_ids > 0) {
        stop_msg <- paste0(
            "No SNPs (inferred as RSIDs) in the dataset start with",
            "'rs'. If these IDs are just some arbitrary value ",
            "rather than RSIDs, set `snp_ids_are_rs_ids=FALSE`"
        )
        if (num_bad_ids == nrow(sumstats_dt)) {
            stop(stop_msg)
        }
        msg <- paste0(
            formatC(num_bad_ids, big.mark = ","), " SNPs (inferred as ",
            "RSIDs) don't start with 'rs' and will be removed"
        )
        message(msg)
        # If user wants log, save it to there
        if (log_folder_ind) {
            name <- "snp_missing_rs"
            name <- get_unique_name_log_file(name = name,
                                             log_files = log_files)
            write_sumstats(
                sumstats_dt = sumstats_dt[!grep("^rs", SNP), ],
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
        sumstats_dt <- sumstats_dt[grep("^rs", SNP), ]

        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}
