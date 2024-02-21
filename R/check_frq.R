#' Ensure all SNPs have frq score above threshold
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
check_frq <- function(sumstats_dt, path, FRQ_filter, log_folder_ind,
                      check_save_out, tabix_index, nThread, log_files) {
    FRQ <- NULL
    col_headers <- names(sumstats_dt)
    if ("FRQ" %in% col_headers && FRQ_filter > 0) {
        message("Filtering SNPs based on FRQ.")
        # use data table for speed
        num_bad_ids <- nrow(sumstats_dt[FRQ < FRQ_filter, ])
        if (num_bad_ids > 0) {
            msg <- paste0(
                formatC(num_bad_ids, big.mark = ","), " SNPs",
                " are below the FRQ threshold of ", FRQ_filter,
                " and will be removed."
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "frq_filter"
                name <- get_unique_name_log_file(
                    name = name,
                    log_files = log_files
                )
                write_sumstats(
                    sumstats_dt = sumstats_dt[FRQ < FRQ_filter, ],
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
                    paste0(
                        check_save_out$log_folder, "/", name,
                        check_save_out$extension
                    )
            }
            sumstats_dt <- sumstats_dt[FRQ >= FRQ_filter, ]
        }
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}
