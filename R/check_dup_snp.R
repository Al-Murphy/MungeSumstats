#' Ensure all rows have unique SNP IDs, drop those that don't
#'
#' @inheritParams format_sumstats  
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data 
#' table object and log files list
#' @keywords internal
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_snp <- function(sumstats_dt, path,log_folder_ind, check_save_out,
                          tabix_index, nThread, log_files){
  SNP = NULL
  col_headers <- names(sumstats_dt)
  if("SNP" %in% col_headers){
    message("Checking for duplicate SNPs from SNP ID.")
    # Try to remove duplicated RSIDs
    data.table::setkey(sumstats_dt,SNP)
    dups <- duplicated(sumstats_dt, by = data.table::key(sumstats_dt))
    if(sum(dups)>0){
      msg <- paste0(formatC(sum(dups),big.mark = ",")," RS IDs are duplicated ",
                    "in the sumstats file. These duplicates will be removed")
      message(msg)
      #If user wants log, save it to there
      if(log_folder_ind){
        name <- "dup_snp_id"
        name <- get_unique_name_log_file(name=name,log_files=log_files)
        write_sumstats(sumstats_dt = sumstats_dt[dups,],
                       save_path=
                         paste0(check_save_out$log_folder,
                                "/",name,
                                check_save_out$extension),
                       sep=check_save_out$sep,
                       tabix_index = tabix_index,
                       nThread = nThread)
        log_files[[name]] <- 
          paste0(check_save_out$log_folder,"/",name,check_save_out$extension)
      } 
      sumstats_dt <- unique(sumstats_dt, by = data.table::key(sumstats_dt))
  
      return(list("sumstats_dt"=sumstats_dt,"log_files"=log_files))
    }
  }
  return(list("sumstats_dt"=sumstats_dt,"log_files"=log_files))
}
