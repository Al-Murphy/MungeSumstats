#' Ensure all rows are unique based on SNP,CHR,BP,A1,A2, drop those that aren't
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and log files list
#' @keywords internal
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_row <- function(sumstats_dt,
                         path,
                         log_folder_ind,
                         check_save_out,
                         tabix_index,
                         nThread,
                         log_files) {
  BP <- CHR <- A1 <- A2 <- SNP <- NULL
  col_headers <- names(sumstats_dt)
  #just check for duplicates across all columns or those would need to be 
  #unique
  #This check is fine for Indels and non-bi-allelic SNPs
  dup_cols <- col_headers
  if(sum(c("SNP","CHR","BP","A1","A2") %in% col_headers) == 5)
    dup_cols <- c("SNP","CHR","BP","A1","A2")
  message("Checking for duplicated rows.")
  # Try to remove duplicated Positions
  data.table::setkeyv(sumstats_dt, dup_cols)
  dups <- duplicated(sumstats_dt, by = data.table::key(sumstats_dt))
  if (sum(dups) > 0) {
    msg <- paste0(
      formatC(sum(dups), big.mark = ","),
      " sumstat rows are",
      " duplicated. These duplicates will",
      " be removed."
    )
    message(msg)
    # If user wants log, save it to there
    if (log_folder_ind) {
      name <- "dup_row"
      name <- get_unique_name_log_file(
        name = name,
        log_files = log_files
      )
      write_sumstats(
        sumstats_dt = sumstats_dt[dups, ],
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
    sumstats_dt <- unique(sumstats_dt,
                          by = data.table::key(sumstats_dt)
    )
    
    return(list(
      "sumstats_dt" = sumstats_dt,
      "log_files" = log_files
    ))
  }
  else{
    return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
  }
}