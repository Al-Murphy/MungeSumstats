#' Convert to \code{VRanges}
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @return \code{VRanges} object
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom VariantAnnotation makeVRangesFromGRanges
to_VRanges <- function(sumstats_dt){ 
    gr <- GenomicRanges::makeGRangesFromDataFrame(sumstats_dt, 
                                                  keep.extra.columns = TRUE, 
                                                  seqnames.field = "CHR", 
                                                  start.field = "BP", 
                                                  end.field = "BP")
    message("Converting summary statistics to VRanges.")
    gr$dummy <- "GWAS"
    vr <- VariantAnnotation::makeVRangesFromGRanges(gr, 
                                                    ref.field = "A1", alt.field = "A2", 
                                                    keep.extra.columns = TRUE,
                                                    sampleNames.field = "dummy")
    return(vr)
}