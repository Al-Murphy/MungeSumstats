#' Standardize the CHR column
#'
#' Maps chromosome names to the default Ensembl/NCBI naming style and removes
#' SNPs with nonstandard CHR entries. Optionally, also removes SNPs on
#' user-specified chromosomes.
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
                      nThread,
                      tabix_index,
                      log_folder_ind) {
  CHR <- NULL

  # The CHR column needs to be a character vector for gsub substitution to work
  sumstats_dt[, CHR := as.character(CHR)]

  ### Reformat chromosome names according to the default style (Ensembl/NCBI)
  # Remove the "chr" prefix
  sumstats_dt[, CHR := gsub("chr", "", CHR, ignore.case = TRUE)]
  # Remove the "ch" prefix
  sumstats_dt[, CHR := gsub("ch", "", CHR, ignore.case = TRUE)]
  # Rename "23" to "X"
  sumstats_dt[, CHR := gsub("23", "X", CHR)]
  # Rename "M" to "MT"
  sumstats_dt[, CHR := gsub("M", "MT", CHR, ignore.case = TRUE)]
  # Make all chromosome names uppercase
  sumstats_dt[, CHR := toupper(CHR)]

  ### Remove rows with nonstandard CHR entries
  standard_chrs <- c(1:22, "X", "Y", "MT")
  nonstandard_rows <- which(!(sumstats_dt$CHR %in% standard_chrs))
  if (length(nonstandard_rows) > 0L) {
    message("Removing ",
            formatC(length(nonstandard_rows), big.mark = ","),
            " SNPs with nonstandard CHR entries.")
  }

  ### Remove SNPs on user-specified chromosomes, if requested
  rmv_chr_rows <- c()
  if (!is.null(rmv_chr)) {
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
    save_path <- paste0(check_save_out$log_folder,
                        "/",
                        name,
                        check_save_out$extension)

    write_sumstats(
      sumstats_dt = sumstats_dt[all_removed_rows],
      save_path = save_path,
      sep = check_save_out$sep,
      #don't tab indx as could be miss values & cause err
      #tabix_index = tabix_index,
      nThread = nThread
    )
    log_files[[name]] <- save_path
  }

  # Remove the SNPs identified above, if any
  sumstats_dt <- sumstats_dt[!all_removed_rows]

  return(list(sumstats_dt = sumstats_dt, log_files = log_files))
}
