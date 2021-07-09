#' Convert to \code{GRanges}
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @return \code{GRanges} object
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
to_GRanges <- function(sumstats_dt){
    message("Converting summary statistics to Genomic Ranges.")
    gr <- GenomicRanges::makeGRangesFromDataFrame(sumstats_dt, 
                                                  keep.extra.columns = TRUE, 
                                                  seqnames.field = "CHR", 
                                                  start.field = "BP", 
                                                  end.field = "BP")
    return(gr)
}