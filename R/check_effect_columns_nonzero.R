#' Ensure that the standard error (se) is positive for all SNPs
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
check_effect_columns_nonzero <- function(sumstats_dt, path,
                                         effect_columns_nonzero,
                                         log_folder_ind, check_save_out,
                                         tabix_index, nThread, log_files) {
    .SD <- NULL
    col_headers <- names(sumstats_dt)
    effect_columns <- c("BETA", "OR", "LOG_ODDS", "SIGNED_SUMSTAT")
    if (sum(effect_columns %in% col_headers) >= 1 && effect_columns_nonzero) {
        message("Filtering effect columns, ensuring none equal 0.")
        # filter to effect columns in the data
        effect_columns_dat <- effect_columns[effect_columns %in% col_headers]
        # ensure numeric
        sumstats_dt[, (effect_columns_dat) := lapply(.SD, as.numeric),
            .SDcols = effect_columns_dat
        ]
        # check if any equal 0 - use data table for speed
        bad_ids <- sumstats_dt[, Reduce(`|`, lapply(.SD, `==`, 0)),
            .SDcols = effect_columns_dat
        ]
        num_bad_ids <- sum(bad_ids)
        if (num_bad_ids > 0) {
            msg <- paste0(
                formatC(num_bad_ids, big.mark = ","), " SNPs",
                " have effect values = 0 and will be removed"
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "effect_col_zero"
                name <- get_unique_name_log_file(name = name, log_files = log_files)
                write_sumstats(
                    sumstats_dt = sumstats_dt[bad_ids, ],
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
            sumstats_dt <- sumstats_dt[!bad_ids]
        }
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}
