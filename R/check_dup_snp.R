#' Ensure all rows have unique SNP IDs, drop those that don't
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and log files list
#' @keywords internal
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_snp <- function(sumstats_dt,
                          indels,
                          path,
                          log_folder_ind,
                          check_save_out,
                          tabix_index,
                          nThread,
                          log_files,
                          bi_allelic_filter,
                          check_dups) {
    SNP <- NULL
    col_headers <- names(sumstats_dt)
    #only remove dups if bi-allelic filter selected, RS IDs are not unique for
    #non bi-allelic SNPs
    if ("SNP" %in% col_headers && bi_allelic_filter && check_dups) {
        message("Checking for duplicate SNPs from SNP ID.")
        #remove indels from this check if sleected since these can have same pos
        indel_dt <- data.table::data.table()
        if(sum(c("A1", "A2") %in% col_headers) == 2 & indels){
          #identify Indels based on num char in A1, A2
          num_indels <- nrow(sumstats_dt[(nchar(A1)>1 | nchar(A2)>1),])
          if(num_indels>0){
            msg <- paste0("Found ",
                          formatC(num_indels,big.mark = ","),
                          " Indels. These won't",
                          " be checked for duplicates based on RS ID ",
                          "as there can be multiples.",
                          "\nWARNING If your sumstat ",
                          "doesn't contain Indels, set the ",
                          "indel param to FALSE & rerun ",
                          "MungeSumstats::format_sumstats()")
            message(msg)
            indel_dt <- sumstats_dt[(nchar(A1)>1 | nchar(A2)>1),]
            #update sumstats to excl indels for check, add back later
            sumstats_dt <- sumstats_dt[!(nchar(A1)>1 | nchar(A2)>1),]
          }
        }
        # Try to remove duplicated RSIDs
        data.table::setkey(sumstats_dt, SNP)
        dups <- duplicated(sumstats_dt,
            by = data.table::key(sumstats_dt)
        )
        if (sum(dups) > 0) {
            msg <- paste0(
                formatC(sum(dups), big.mark = ","), " RSIDs are duplicated ",
                "in the sumstats file. These duplicates will be removed"
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "dup_snp_id"
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
            #join back indels
            sumstats_dt <-
              data.table::rbindlist(
                list(sumstats_dt, indel_dt))
            return(list(
                "sumstats_dt" = sumstats_dt,
                "log_files" = log_files
            ))
        }
        #join back indels
        sumstats_dt <-
          data.table::rbindlist(
            list(sumstats_dt, indel_dt))
    }
    return(list(
        "sumstats_dt" = sumstats_dt,
        "log_files" = log_files
    ))
}
