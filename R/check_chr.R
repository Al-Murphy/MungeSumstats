#' Standardize the CHR column
#'
#' Renames "23" to "X", makes X/Y/MT uppercase and the "chr" prefix lowercase,
#' and removes SNPs with nonstandard CHR entries. Optionally, also removes the
#' "chr" prefix and SNPs on user-specified chromosomes.
#'
#' @param sumstats_dt data.table with summary statistics
#' @param log_files list of locations for all log files
#' @param check_save_out list of parameters for saved files
#' @inheritParams format_sumstats
#' @return list containing the updated summary statistics data.table and the
#'   updated log file locations list
#' @keywords internal
check_chr <- function(sumstats_dt,
                      log_files,
                      check_save_out,
                      rmv_chr,
                      rmv_chrPrefix,
                      nThread,
                      tabix_index,
                      log_folder_ind) {
  CHR <- NULL

  ### Sometimes X is labeled as 23
  sumstats_dt[, CHR := gsub("23", "X", CHR)]

  ### Make X/Y/MT uppercase
  sumstats_dt[, CHR := gsub("x", "X", CHR)]
  sumstats_dt[, CHR := gsub("y", "Y", CHR)]
  sumstats_dt[, CHR := gsub("mt", "MT", CHR)]

  ### If specified, remove the "chr" prefix
  if (rmv_chrPrefix) {
    message("Removing 'chr' prefix from CHR.")
    sumstats_dt[, CHR := gsub("chr", "", CHR, ignore.case = TRUE)]
    standard_chrs <- c(1:22, "X", "Y", "MT")
  } else {
    ### Otherwise, make the "chr" prefix lowercase
    sumstats_dt[, CHR := gsub("CHR", "chr", CHR)]
    standard_chrs <- c(paste0("chr", 1:22), "X", "Y", "MT")
  }

  ### Remove rows with nonstandard CHR entries
  nonstandard_rows <- which(!(sumstats_dt$CHR %in% standard_chrs))
  if (length(nonstandard_rows) > 0L) {
    message(
      "Removing ",
      formatC(length(nonstandard_rows), big.mark = ","),
      " SNPs with nonstandard CHR entries."
    )
  }

  ### If specified, remove SNPs on specific chromosomes
  rmv_chr_rows <- c()
  if (!is.null(rmv_chr)) {
    # Standardize user-specified chromosomes
    rmv_chr <- gsub("23", "X", rmv_chr)
    rmv_chr <- gsub("x", "X", rmv_chr)
    rmv_chr <- gsub("y", "Y", rmv_chr)
    rmv_chr <- gsub("mt", "MT", rmv_chr)
    if (rmv_chrPrefix) {
      rmv_chr <- gsub("chr", "", rmv_chr, ignore.case = TRUE)
    } else {
      rmv_chr <- gsub("CHR", "chr", rmv_chr)
    }

    # Check for chromosomes to be removed
    rmv_chr_rows <- which(sumstats_dt$CHR %in% rmv_chr)
    if (length(rmv_chr_rows) > 0L) {
      message(
        formatC(length(rmv_chr_rows), big.mark = ","),
        " SNPs are on chromosomes ",
        paste(rmv_chr, collapse = ", "),
        " and will be removed."
      )
    }
  }

  # Vector of row numbers for all removed SNPs
  all_removed_rows <- sort(unique(c(nonstandard_rows, rmv_chr_rows)))

  ### Save a log of removed SNPs if the user wants it
  if (log_folder_ind && (length(all_removed_rows) > 0L)) {
    name <- "chr_excl"
    name <- get_unique_name_log_file(name = name,
                                     log_files = log_files)
    save_path <- paste0(
      check_save_out$log_folder,
      "/",
      name,
      check_save_out$extension
    )

    write_sumstats(sumstats_dt = sumstats_dt[all_removed_rows],
                   save_path = save_path,
                   sep = check_save_out$sep,
                   tabix_index = tabix_index,
                   nThread = nThread)
    log_files[[name]] <- save_path
  }

  # Remove the SNPs identified above, if any
  sumstats_dt <- sumstats_dt[!all_removed_rows]

  return(list(sumstats_dt = sumstats_dt, log_files = log_files))
}
