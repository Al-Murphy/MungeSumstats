#' GenomicRanges to data.table
#' 
#' Convert a \link[GenomicRanges]{GRanges} into a \link[data.table]{data.table}.
#' @param gr A \link[GenomicRanges]{GRanges} object.
#' 
#' @return A data.table object.
#' @source 
#' \href{https://rdrr.io/github/databio/GenomicDistributions/src/R/utility.R}{
#' Code adapted from GenomicDistributions.}
#' @keywords internal
#' @importFrom GenomicRanges seqnames start end elementMetadata
#' @importFrom data.table data.table
granges_to_dt  <- function(gr) {
    if(is.null(gr)) return(gr)
    #### Convert metadata ####
    DF <- GenomicRanges::elementMetadata(gr)
    #### Combine metadata with ranges data ####
    if( ncol(DF) > 0) { 
        meta <- DF_to_dt(DF = DF)
        DT <- data.table::data.table(
            chr=as.vector(GenomicRanges::seqnames(gr)), 
            start=GenomicRanges::start(gr),
            end=GenomicRanges::end(gr), 
            meta)
    } else {
        DT <- data.table::data.table(
            chr=as.vector(GenomicRanges::seqnames(gr)),
            start=GenomicRanges::start(gr),
            end=GenomicRanges::end(gr))
    }
    DT <- check_numeric(DT)
    return(DT)
}