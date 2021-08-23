#' Ensure that the standard error (se) is positive for all SNPs
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data 
#' table object and the log file list
#' @keywords internal
check_pos_se <- function(sumstats_dt, path, pos_se,log_folder_ind, 
                         check_save_out, tabix_index, nThread, log_files){
    SE = .SD = NULL
    col_headers <- names(sumstats_dt)
    if("SE" %in% col_headers && pos_se){
        message("Filtering SNPs, ensuring SE>0.")
        #use data table for speed
        num_bad_ids <- nrow(sumstats_dt[SE<=0,])
        if(num_bad_ids>0){
            msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs",
                          " have SE values <= 0 and will be removed")
            message(msg)
            #If user wants log, save it to there
            if(log_folder_ind){
                name <- "se_neg"
                name <- get_unique_name_log_file(name=name,log_files=log_files)
                write_sumstats(sumstats_dt = sumstats_dt[SE<=0,],
                               save_path=
                                   paste0(check_save_out$log_folder,
                                          "/",name,
                                          check_save_out$extension),
                               sep=check_save_out$sep,
                               tabix_index = tabix_index,
                               nThread = nThread)
                log_files[[name]] <- 
                    paste0(check_save_out$log_folder,"/",name,
                            check_save_out$extension)
            } 
            sumstats_dt <- sumstats_dt[SE>0,]
        }
        return(list("sumstats_dt"=sumstats_dt,"log_files"=log_files))
    }
    else{
        return(list("sumstats_dt"=sumstats_dt,"log_files"=log_files))
    }
}