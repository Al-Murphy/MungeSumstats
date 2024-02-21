#' Ensure all SNPs have N less than X std dev below mean
#'
#' In case some SNPs were genotyped by a specialized genotyping array and
#' have substantially more samples than others. These will be removed.
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
#' @importFrom stats sd
check_n_num <- function(sumstats_dt, 
                        path, 
                        N_std, 
                        N_dropNA = FALSE,
                        log_folder_ind,
                        check_save_out, 
                        tabix_index, 
                        nThread, 
                        log_files) {
    message("Ensuring all SNPs have N<", N_std, " std dev above mean.")
    N <- NULL
    col_headers <- names(sumstats_dt)
    if ("N" %in% col_headers && N_std > 0) {
        mean_N <- mean(sumstats_dt$N)
        sd_N <- stats::sd(sumstats_dt$N)
        num_bad_ids <- nrow(sumstats_dt[N > ((N_std * sd_N) + mean_N), ])
        if (num_bad_ids > 0) {
            msg <- paste0(
                formatC(num_bad_ids, big.mark = ","), " SNPs have N values ",
                N_std, " standard deviations above the mean",
                " and will be removed"
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "n_large"
                name <- get_unique_name_log_file(name = name,
                                                 log_files = log_files)
                write_sumstats(
                    sumstats_dt = sumstats_dt[N > ((N_std * sd_N) + mean_N), ],
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
            sumstats_dt <- sumstats_dt[N <= ((N_std * sd_N) + mean_N), ]
        }
        if (!N_dropNA) {
            message("Removing rows where is.na(N)")
            n_NAs <- sum(is.na(sumstats_dt$N))
            message(
                formatC(n_NAs, big.mark = ","),
                " SNPs have N values that are NA and will be removed."
            )
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "n_null"
                name <- get_unique_name_log_file(name = name,
                                                 log_files = log_files)
                write_sumstats(
                    sumstats_dt = sumstats_dt[!complete.cases(N)],
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
            sumstats_dt <- sumstats_dt[complete.cases(N)]
        }
        return(list("sumstats_dt" = sumstats_dt,
                    "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt,
                    "log_files" = log_files))
    }
}
