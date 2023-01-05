#' To \code{GRanges}
#' 
#' 
#' Convert a \link[data.table]{data.table} to \link[GenomicRanges]{GRanges}.
#' @param sumstats_dt data table obj of the summary statistics file 
#' for the GWAS.
#' @param style \code{GRanges} style to convert to, "NCBI" or "UCSC".
#' @inheritParams GenomicRanges::makeGRangesFromDataFrame
#' @return \code{GRanges} object
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevelsStyle
to_granges <- function(sumstats_dt,
                       seqnames.field = "CHR",
                       start.field = "BP",
                       end.field = "BP",
                       style = c("NCBI", "UCSC")) {
    message("Converting summary statistics to GenomicRanges.")
    gr <- GenomicRanges::makeGRangesFromDataFrame(
        df = sumstats_dt,
        keep.extra.columns = TRUE,
        seqnames.field = seqnames.field,
        start.field = start.field,
        end.field = end.field
    )
    
    suppressMessages(suppressWarnings(
        GenomeInfoDb::seqlevelsStyle(gr) <- style[1]
    ))
    return(gr)
}
