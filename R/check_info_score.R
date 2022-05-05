#' Ensure all SNPs have info score above threshold
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations. 
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
check_info_score <- function(sumstats_dt,
                             INFO_filter, 
                             log_folder_ind,
                             check_save_out, 
                             tabix_index, 
                             nThread, 
                             log_files) {
    INFO <- NULL
    col_headers <- names(sumstats_dt)
    if ("INFO" %in% col_headers && INFO_filter > 0) {
        msg1 <- "Filtering SNPs based on INFO score."
        message(msg1)
        # use data table for speed
        num_bad_ids <- nrow(sumstats_dt[INFO < INFO_filter, ])
        if (num_bad_ids > 0) {
            msg2 <- paste0(
                formatC(num_bad_ids, big.mark = ","), " SNPs",
                " are below the INFO threshold of ", INFO_filter,
                " and will be removed."
            )
            message(msg2)
            #### If user wants log, save it to there ####
            if (log_folder_ind) {
                log_out <- check_info_score_log(sumstats_dt=sumstats_dt,
                                                log_files=log_files,
                                                INFO_filter=INFO_filter,
                                                check_save_out=check_save_out,
                                                tabix_index=tabix_index,
                                                nThread=nThread)
                sumstats_dt <- log_out$sumstats_dt
                log_files <- log_out$log_files
            }
            sumstats_dt <- sumstats_dt[INFO >= INFO_filter, ]
        } else {
            msg2b <- paste0("All rows have INFO>=",INFO_filter)
            message(msg2b)
        }
    } else {
        if(!"INFO" %in% col_headers){
            msg3 <- "INFO column not available. Skipping INFO score filtering step."
            message(msg3)
        } else if (INFO_filter == 0){
            msg3 <- "INFO_filter==0. Skipping INFO score filtering step."
            message(msg3)
        }  
    }
    return(list("sumstats_dt" = sumstats_dt,
                "log_files" = log_files))
}