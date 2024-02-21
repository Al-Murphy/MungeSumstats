#' Drop Indels from summary statistics
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @inheritParams format_sumstats
#' @source 
#' \code{
#' sumstats_dt <- MungeSumstats:::formatted_example()
#' sumstats <- check_drop_indels(sumstats_dt = sumstats_dt, 
#'                               drop_indels = TRUE)

#' }
#' @returns list containing sumstats_dt, 
#' the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table :=
check_drop_indels <- function(sumstats_dt, 
                              drop_indels,
                              path,
                              log_folder_ind,
                              check_save_out,
                              tabix_index,
                              nThread,
                              log_files) {
    A1 <- A2 <- NULL
    # Make sure A1 and A2 are present
    col_headers <- names(sumstats_dt)
    if (sum(c("A1", "A2") %in% col_headers) == 2 & drop_indels) {
        # Check for indels
        #identify Indels based on num char in A1, A2
        indels_dt <- sumstats_dt[(nchar(A1)>1 | nchar(A2)>1),]
        num_indels <- nrow(indels_dt)
        if (num_indels > 0) {
            msg <- paste0("Found ",
                        formatC(num_indels,big.mark = ","),
                        " Indels. These will",
                        " be removed from the sumstats. ",
                        "\nWARNING If you want ",
                        "to keep Indels, set the ",
                        "drop_indel param to FALSE & rerun ",
                        "MungeSumstats::format_sumstats()")
            message(msg)
            #remove the indels
            sumstats_dt <- sumstats_dt[!(nchar(A1)>1 | nchar(A2)>1),]
            # If user wants log, save it to there
            if (log_folder_ind) {
              name <- "indel"
              name <- get_unique_name_log_file(
                name = name,
                log_files = log_files
              )
              write_sumstats(
                sumstats_dt = indels_dt,
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
        }
    }    
    return(list("sumstats_dt" = sumstats_dt,
                "log_files" = log_files))
}
