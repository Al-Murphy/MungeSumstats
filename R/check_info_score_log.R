check_info_score_log <- function(sumstats_dt,
                                 log_files,
                                 INFO_filter,
                                 check_save_out,
                                 tabix_index,
                                 nThread){
    INFO <- NULL
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
    return(list(log_files=log_files,
                sumstats_dt=sumstats_dt))
}
