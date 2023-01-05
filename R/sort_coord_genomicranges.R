#' Sort sum stats: GenomicRanges
#' 
#' Sort summary statistics table by genomic coordinates using a slower 
#' (but in some cases more robust) \code{GenomicRanges} strategy 
#' @inheritParams sort_coords
#' @returns Sorted sumstats_dt
#' 
#' @keywords internal
#' @importFrom GenomicRanges sort.GenomicRanges
#' @importFrom data.table setnames :=
sort_coord_genomicranges <- function(sumstats_dt){
    end <- NULL;
    ### sort.GenomicRanges is less efficient than setorderv,
    ### but more consistently reliable for tabix-indexing.
    keep_end <- "end" %in% tolower(names(sumstats_dt))
    gr <- to_granges(sumstats_dt)
    gr <- GenomicRanges::sort.GenomicRanges(gr)
    sumstats_dt <- granges_to_dt(gr)
    remove(gr)
    if(isFALSE(keep_end)) sumstats_dt[,end:=NULL]
    data.table::setnames(sumstats_dt, c("chr","start"), c("CHR","BP"))
    return(sumstats_dt)
}