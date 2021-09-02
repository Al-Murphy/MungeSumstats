#' Convert data.table to GRanges object
#'
#' @param dat data.table of genomic data
#' @param style conversion to GRanges style, NCBI or UCSC
#' @param chrom_col column for chromosome
#' @param start_col column for start position
#' @param end_col column for end position (is the same as start for snps)
#' @return GRanges object of the genomic data
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevelsStyle
dt_to_granges <- function(dat,
                          style = c("NCBI", "UCSC"),
                          chrom_col = "CHR",
                          start_col = "POS",
                          end_col = "POS") {
    dat[["SEQnames"]] <- dat[[chrom_col]]
    gr.snp <- 
        GenomicRanges::makeGRangesFromDataFrame(dat,
                                                seqnames.field = "SEQnames",
                                                start.field = start_col,
                                                end.field = end_col,
                                                keep.extra.columns = TRUE)
    suppressMessages(suppressWarnings(
        GenomeInfoDb::seqlevelsStyle(gr.snp) <- style[1]
    ))
    return(gr.snp)
}