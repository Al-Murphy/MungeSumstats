#' Ensure that the standard error (se) is positive for all SNPs
#' Also impute se if missing
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
check_pos_se <- function(sumstats_dt, path, pos_se, log_folder_ind,
                         check_save_out, tabix_index, nThread, log_files,
                         impute_se) {
    `%nin%` = negate(`%in%`)
    SE <- .SD <- NULL
    col_headers <- names(sumstats_dt)
    no_imp_msg <- "SE is not present but can be imputed with "
    no_imp_msg2 <- ". Set impute_se=TRUE and rerun to do this."
    #message about the last way to impute beta from other columns
    imp_cols <- FALSE
    se_imputed <- FALSE
    
    if ("SE" %nin% col_headers) {
        
        derive_msg = message("The sumstats SE column is not present...")
        
        if ("Z" %in% col_headers & "BETA" %in% col_headers) {
            imp_cols <- "Z & BETA"
            if(impute_se){
                message(paste0(derive_msg,"Deriving SE from Z and BETA"))
                sumstats_dt[,SE := BETA / Z]
                se_imputed <- TRUE
            }  
        } else if ("BETA" %in% col_headers & "P" %in% col_headers) {
            imp_cols <- "BETA & P"
            if(impute_se){
                # https://www.biostars.org/p/431875/
                message(paste0(derive_msg,"Deriving SE from Beta and P"))
                sumstats_dt[,SE := abs(BETA/ qnorm(P/2))]
                se_imputed <- TRUE
            }  
        }
        #tell the user if they could impute beta but didn't because of input param
        if(!impute_se && !isFALSE(imp_cols)){
            no_imp_msg <- paste0(no_imp_msg,imp_cols,no_imp_msg2)
            message(no_imp_msg)
        }
        
        # If user wants log, save it to there
        if (log_folder_ind && nrow(sumstats_dt[is.na(SE), ])>0) {
            name <- "se_na"
            name <- get_unique_name_log_file(name = name,
                                             log_files = log_files)
            write_sumstats(
                sumstats_dt =
                    sumstats_dt[is.na(SE),],
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
                paste0(check_save_out$log_folder, "/",
                       name, check_save_out$extension)
        }
        
        # if user specifies add a column to notify of the imputation
        if (imputation_ind && se_imputed) {
            sumstats_dt[, IMPUTATION_SE := TRUE]
        }
        
        ## Remove NA values if introduced
        if (impute_se | "SE" %in% col_headers ){
            sumstats_dt <- sumstats_dt[!is.na(SE), ]
        } 
        
    }
    
    if ("SE" %in% col_headers && pos_se) {
        message("Filtering SNPs, ensuring SE>0.")
        # use data table for speed
        num_bad_ids <- nrow(sumstats_dt[SE <= 0, ])
        if (num_bad_ids > 0) {
            msg <- paste0(
                formatC(num_bad_ids, big.mark = ","), " SNPs",
                " have SE values <= 0 and will be removed"
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "se_neg"
                name <- get_unique_name_log_file(name = name,
                                                 log_files = log_files)
                write_sumstats(
                    sumstats_dt = sumstats_dt[SE <= 0, ],
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
            sumstats_dt <- sumstats_dt[SE > 0, ]
        }
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}