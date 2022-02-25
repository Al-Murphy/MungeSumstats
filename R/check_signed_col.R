#' Ensure that there is at least one signed column in summary statistics file
#' Impute beta if user requests
#'
#' @param sumstats_dt data table obj of the summary statistics
#' file for the GWAS
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return null
#' @keywords internal

check_signed_col <-
    function(sumstats_dt,impute_beta, log_folder_ind, rsids, imputation_ind,
             check_save_out,log_files,nThread) {
        col_headers <- names(sumstats_dt)
        signed_stat_column_names <- c("Z", "OR", "BETA",
                                      "LOG_ODDS", "SIGNED_SUMSTAT")
        `%nin%` = negate(`%in%`)
        stp_msg <- paste0(
            "ERROR: cannot find a column name representing signed ",
            "statistic in GWAS sumstats file:\n",
            "'Z','OR','BETA','LOG_ODDS','SIGNED_SUMSTAT'"
        )
        
        no_imp_msg <- "BETA is not present but can be imputed with "
        no_imp_msg2 <- ". Set impute_beta=TRUE and rerun to do this."
        msg <- "The sumstats Beta column is not present... "
        #message about the last way to impute beta from other columns
        imp_cols <- FALSE
        beta_imputed <- FALSE
        #impute BETA from other columns
        if ("BETA" %nin% col_headers) {
            if ("OR" %in% col_headers) {
                imp_cols <- "OR"
                if(impute_beta){
                    message(paste0(msg,"Deriving BETA from OR"))
                    sumstats_dt[,BETA := log(OR)]
                    beta_imputed <- TRUE
                }
            } else if ("Z" %in% col_headers & "SE" %in% col_headers) {
                imp_cols <- "Z & SE"
                if(impute_beta){
                    message(paste0(msg,"Deriving BETA from Z and SE"))
                    sumstats_dt[,BETA := Z * SE]
                    beta_imputed <- TRUE
                }  
            } else if ("Z" %in% col_headers && "N" %in% col_headers &&
                       "FRQ" %in% col_headers){
                imp_cols <- "Z, N & FRQ"
                # https://www.biostars.org/p/319584/
                # https://www.nature.com/articles/ng.3538
                if(impute_beta){
                    message(paste0(msg,"Deriving BETA from Z, N, and FRQ"))
                    sumstats_dt[,BETA := Z/sqrt(2*FRQ*(1-FRQ)*(N+Z^2))]
                    beta_imputed <- TRUE
                }  
            } else if ("Z" %in% col_headers && "N" %in% col_headers &&
                       "P" %in% col_headers){
                imp_cols <- "Z, N & P"
                if(impute_beta){
                    message(paste0(msg,"Deriving BETA from Z, N, and P"))
                    sumstats_dt[,BETA := Z/sqrt(qchisq(P, N))]
                    beta_imputed <- TRUE
                }  
            } else if (sum(signed_stat_column_names %in% col_headers) < 1) {
                stop(stp_msg)
            }
            #tell the user if they could impute beta but didn't because of input param
            if(!impute_beta && !isFALSE(imp_cols)){
                no_imp_msg <- paste0(no_imp_msg,imp_cols,no_imp_msg2)
                message(no_imp_msg)
            }
            
            # If user wants log, save it to there
            if (log_folder_ind && nrow(sumstats_dt[is.na(BETA), ])>0) {
                name <- "beta_na"
                name <- get_unique_name_log_file(name = name,
                                                 log_files = log_files)
                write_sumstats(
                    sumstats_dt =
                        sumstats_dt[is.na(BETA),],
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
            if (imputation_ind && beta_imputed) {
                sumstats_dt[, IMPUTATION_BETA := TRUE]
            }
            
            ## Remove NA values if introduced
            sumstats_dt <- sumstats_dt[!is.na(BETA), ]
            
            return(list(
                "sumstats_dt" = sumstats_dt,
                "rsids" = rsids, "log_files" = log_files
            ))
        }
        return(list(
            "sumstats_dt" = sumstats_dt,
            "rsids" = rsids, "log_files" = log_files
        ))
    }