#' Remove non-standard VCF columns
#'
#' @return sumstats_file
#'
#' @inheritParams format_sumstats
#' @param keep_extra_cols logical, should extra, non-essential columns from
#' input be kept
#' @importFrom VariantAnnotation scanVcfHeader
#' @importFrom data.table :=
#' @keywords internal
remove_nonstandard_vcf_cols <- function(sample_id,
                                        sumstats_file,
                                        keep_extra_cols = FALSE,
                                        standardise_headers = TRUE) {
    if(standardise_headers){
        sumstats_file <- standardise_sumstats_column_headers_crossplatform(
            sumstats_dt = sumstats_file)$sumstats_dt
    }
    keep_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "INFO", sample_id,
                   sumstatsColHeaders$Corrected)
    colsToRemove <- names(sumstats_file)[!names(sumstats_file) %in% keep_cols]
    if (!is.null(colsToRemove) && keep_extra_cols == FALSE) {
        message(
            "Removing non-standard columns: ",
            paste(colsToRemove, collapse = ", ")
        )
        sumstats_file[, (colsToRemove) := NULL]
    }
    return(sumstats_file)
}
