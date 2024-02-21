#' Ensure all rows have unique positions, drop those that don't
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and log files list
#' @keywords internal
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_bp <- function(sumstats_dt,
                         bi_allelic_filter,
                         check_dups,
                         indels,
                         path,
                         log_folder_ind,
                         check_save_out,
                         tabix_index,
                         nThread,
                         log_files) {
    BP <- CHR <- A1 <- A2 <- SNP <- NULL
    col_headers <- names(sumstats_dt)
    #non-bi-allelic SNPs can share same BP position so don't check if selected
    if (sum(c("BP", "CHR") %in% col_headers) == 2 && 
        bi_allelic_filter && check_dups) {
        message("Checking for SNPs with duplicated base-pair positions.")
        # Try to remove duplicated Positions
        data.table::setkey(sumstats_dt, BP, CHR)
        #remove indels from this check if sleected since these can have same pos
        indel_dt <- data.table::data.table()
        if(sum(c("A1", "A2") %in% col_headers) == 2 & indels ){
            #identify Indels based on num char in A1, A2
            num_indels <- nrow(sumstats_dt[(nchar(A1)>1 | nchar(A2)>1),])
            if(num_indels>0){
                msg <- paste0("Found ",
                              formatC(num_indels,big.mark = ","),
                              " Indels. These won't",
                              " be checked for duplicates based on base-pair ",
                              "position as there can be multiples.",
                              "\nWARNING If your sumstat ",
                              "doesn't contain Indels, set the ",
                              "indel param to FALSE & rerun ",
                              "MungeSumstats::format_sumstats()")
                message(msg)
                indel_dt <- sumstats_dt[(nchar(A1)>1 | nchar(A2)>1),]
                #run generic dup check on indels
                indel_dt_rtn <- check_dup_row(
                    sumstats_dt = indel_dt,
                    check_dups = check_dups,
                    path = path,
                    log_folder_ind = log_folder_ind,
                    check_save_out = check_save_out,
                    tabix_index = tabix_index,
                    nThread = nThread,
                    log_files = log_files
                )
                indel_dt <- indel_dt_rtn$sumstats_dt
                # update values
                log_files <- indel_dt_rtn$log_files
                sumstats_dt <- sumstats_dt[!(nchar(A1)>1 | nchar(A2)>1),]
            }
        }
        dups <- duplicated(sumstats_dt, by = data.table::key(sumstats_dt))
        if (sum(dups) > 0) {
            msg <- paste0(
                formatC(sum(dups), big.mark = ","),
                " base-pair positions are",
                " duplicated in the sumstats file. These duplicates will",
                " be removed."
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "dup_base_pair_position"
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
                    #don't tab indx as could be miss values & cause err
                    #tabix_index = tabix_index,
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
    else{
        #run generic dup check
        sumstats_return <- check_dup_row(
            sumstats_dt = sumstats_dt,
            check_dups = check_dups,
            path = path,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        sumstats_dt <- sumstats_return$sumstats_dt
        log_files <- sumstats_return$log_files
    }
    return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
}
