#' Ensure all SNPs have info score above threshold
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
check_info_score <- function(sumstats_dt, path, INFO_filter, log_folder_ind,
                             check_save_out, tabix_index, nThread, log_files) {
    INFO <- NULL
    col_headers <- names(sumstats_dt)
    if ("INFO" %in% col_headers && INFO_filter > 0) {
        message("Filtering SNPs based on INFO score.")
        # use data table for speed
        num_bad_ids <- nrow(sumstats_dt[INFO < INFO_filter, ])
        if (num_bad_ids > 0) {
            msg <- paste0(
                formatC(num_bad_ids, big.mark = ","), " SNPs",
                " are below the INFO threshold of ", INFO_filter,
                " and will be removed."
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "info_filter"
                name <- get_unique_name_log_file(
                    name = name,
                    log_files = log_files
                )
                write_sumstats(
                    sumstats_dt = sumstats_dt[INFO < INFO_filter, ],
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
            sumstats_dt <- sumstats_dt[INFO >= INFO_filter, ]
        }
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}
