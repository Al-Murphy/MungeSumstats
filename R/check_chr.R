#' Ensure all SNPs on specified chromosomes are removed
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
check_chr <- function(sumstats_dt,
                      path,
                      rmv_chr,
                      log_folder_ind,
                      check_save_out,
                      tabix_index,
                      nThread,
                      log_files,
                      make_uppercase = TRUE,
                      rmv_chrPrefix = TRUE) {
    CHR <- NULL
    # If CHR present and user specified chromosome to have SNPs removed
    col_headers <- names(sumstats_dt)
    if ("CHR" %in% col_headers && !is.null(rmv_chr)) {

        ### Sometimes X is labeled as 23
        sumstats_dt[, CHR := gsub("23", "X", CHR)]

        #### Remove chr prefix uppercase ####
        if (rmv_chrPrefix) {
            message("Removing 'chr' prefix from CHR.")
            sumstats_dt[, CHR := gsub("chr", "", CHR, ignore.case = TRUE)]
            rmv_chr <- gsub("chr", "", rmv_chr, ignore.case = TRUE)
        }
        #### Make all CHR uppercase ####
        if (make_uppercase) {
            message("Making X/Y CHR uppercase.")
            sumstats_dt[, CHR := gsub("x|23", "X", CHR)]
            sumstats_dt[, CHR := gsub("y", "Y", CHR)]
        }

        # check for chromosomes to be removed
        ### Standardise chromosomes specified
        rmv_chr <- toupper(rmv_chr)
        if (any(rmv_chr %in% unique(sumstats_dt$CHR))) {
            num_bad_ids <- nrow(sumstats_dt[CHR %in% rmv_chr, ])
            msg <- paste0(
                formatC(num_bad_ids, big.mark = ","),
                " SNPs are on chromosomes ",
                paste(rmv_chr, collapse = ", "),
                " and will be removed"
            )
            message(msg)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "chr_excl"
                name <- get_unique_name_log_file(
                    name = name,
                    log_files = log_files
                )
                write_sumstats(
                    sumstats_dt = sumstats_dt[CHR %in% (rmv_chr), ],
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
            # remove rows on these chromosomes
            sumstats_dt <- sumstats_dt[!CHR %in% (rmv_chr), ]
        }

        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    } else {
        return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
    }
}
