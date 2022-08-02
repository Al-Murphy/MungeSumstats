#' Remove SNPs with missing data
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and a log file list.
#' @keywords internal
#' @importFrom stats complete.cases
check_miss_data <- function(sumstats_dt, path, log_folder_ind, check_save_out,
                            tabix_index, nThread, log_files) {
    message("Checking for missing data.")
    col_headers <- names(sumstats_dt)
    # use data table for speed
    # check for rows missing data to be excluded
    # don't check imputation columns
    #also don't check cols MSS creates - SNP_INFO
    ignore_cols <- c(
        col_headers[grepl("^IMPUTATION_", col_headers)],
        "flipped"["flipped" %in% col_headers],
        col_headers[grepl("^convert_", col_headers)],
        "SNP_INFO"["SNP_INFO"%in% col_headers]
    )
    incl_cols <- names(sumstats_dt)[!names(sumstats_dt) %in% ignore_cols]
    if (nrow(sumstats_dt[!complete.cases(sumstats_dt[, incl_cols,
        with = FALSE
    ]), ]) > 0) {
        n_missing <- nrow(
            sumstats_dt[!complete.cases(
                sumstats_dt[, incl_cols, with = FALSE]
            ), ]
        )
        print(sumstats_dt[!complete.cases(
          sumstats_dt[, incl_cols, with = FALSE]
        ), ])
        msg <- paste0(
            "WARNING: ",
            formatC(n_missing, big.mark = ","),
            " rows in sumstats file are missing data and will ",
            "be removed."
        )
        message(msg)
        # If user wants log, save it to there
        if (log_folder_ind) {
            name <- "missing_data"
            name <- get_unique_name_log_file(
                name = name,
                log_files = log_files
            )
            write_sumstats(
                sumstats_dt =
                    sumstats_dt[!complete.cases(sumstats_dt[, incl_cols,
                        with = FALSE
                    ])],
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
                    check_save_out$log_folder, "/",
                    name, check_save_out$extension
                )
        }
        sumstats_dt <-
            sumstats_dt[complete.cases(sumstats_dt[, incl_cols, with = FALSE])]
        if (nrow(sumstats_dt) == 0) {
            stop_msg <- paste(
                "All SNPs have been filtered out of",
                " your summary statistics dataset"
            )
            stop(stop_msg)
        }

        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}
