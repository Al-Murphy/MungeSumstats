#' Convert to \code{VRanges}
#'
#' @param sumstats_dt data table obj of the summary statistics 
#' file for the GWAS.
#' @return \code{VRanges} object
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom VariantAnnotation makeVRangesFromGRanges
to_vranges <- function(sumstats_dt) {
    gr <- to_granges(sumstats_dt)
    message("Converting summary statistics to VRanges.")
    gr$dummy <- "GWAS"
    vr <- VariantAnnotation::makeVRangesFromGRanges(gr,
        ref.field = "A1", alt.field = "A2",
        keep.extra.columns = TRUE,
        sampleNames.field = "dummy"
    )
    return(vr)
}
