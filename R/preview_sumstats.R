#' Preview formatted sum stats saved to disk
#'
#' Prints the first \code{n} lines of the sum stats.
#'
#' @return No return
#'
#' @inheritParams format_sumstats
#' @keywords internal
#' @importFrom utils capture.output
preview_sumstats <- function(save_path,
                             nrows = 5L) {
    message("Successfully finished preparing sumstats file, preview:")
    vcf_suffixes <- supported_suffixes(tabular = FALSE,
                                       tabular_compressed = FALSE)
    preview <- read_header(
        path = save_path,
        n = nrows,
        skip_vcf_metadata = TRUE)
    message(paste0(capture.output(preview), collapse = "\n"))
}
