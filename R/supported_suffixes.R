#' List supported file formats
#'
#' @param tabular Include tabular formats.
#' @param tabular_compressed Include compressed tabular formats.
#' @param vcf Include Variant Call Format.
#' @param vcf_compressed Include compressed Variant Call Format.
#' @return File formats
#' @keywords internal
supported_suffixes <- function(tabular = TRUE,
                               tabular_compressed = TRUE,
                               vcf = TRUE,
                               vcf_compressed = TRUE) {
    supported <- c()
    suffixes <- c(".tsv", ".txt", ".csv")
    suffixes.gz <- c(paste0(suffixes, ".gz"), paste0(suffixes, ".bgz"))
    suffixes.vcf <- c(".vcf")
    suffixes.vcf.gz <- c(paste0(suffixes.vcf, ".gz"),
                         paste0(suffixes.vcf, ".bgz"))
    if (tabular) supported <- c(supported, suffixes)
    if (tabular_compressed) supported <- c(supported, suffixes.gz)
    if (vcf) supported <- c(supported, suffixes.vcf)
    if (vcf_compressed) supported <- c(supported, suffixes.vcf.gz)
    return(supported)
}
