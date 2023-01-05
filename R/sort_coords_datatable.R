#' Sort sum stats: data.table
#' 
#' Sort summary statistics table by genomic coordinates using a fast
#'  \code{data.table}-native strategy
#' @param chr_col Chromosome column name.
#' @param start_col Genomic start position column name.
#' @param start_col Genomic end position column name.
#' @inheritParams sort_coords
#' @returns Sorted sumstats_dt
#' 
#' @keywords internal
#' @importFrom data.table := setorderv
#' @importFrom GenomeInfoDb rankSeqlevels
sort_coords_datatable <- function(sumstats_dt,
                                  chr_col="CHR",
                                  start_col="BP",
                                  end_col=start_col){
    CHR <- NULL;
    
    #### Check args ####
    if(!chr_col %in% names(sumstats_dt)){
        stp <- "chr_col must be in colnames of sumstats_dt."
        stop(stp)
    }
    if(!start_col %in% names(sumstats_dt)){
        stp <- "start_col must be in colnames of sumstats_dt."
        stop(stp)
    }
    if(!end_col %in% names(sumstats_dt)){
        stp <- "end_col must be in colnames of sumstats_dt."
        stop(stp)
    }
    #### Order seqnames #### 
    ans_seqnames <- sumstats_dt[,chr_col,with=FALSE][[1]]
    seqlevels <- levels(ans_seqnames)
    if (is.null(seqlevels)) {
        seqlevels <- unique(ans_seqnames)
        if (!is.character(seqlevels)) seqlevels <- as.character(seqlevels)
    }
    seqlevels[GenomeInfoDb::rankSeqlevels(seqnames = seqlevels)] <- seqlevels 
    ### Turn CHR into an ordered factor to account for X/Y/MT chroms
    sumstats_dt[, CHR:=factor(CHR, levels = seqlevels, ordered = TRUE)]
    ## IMPORTANT: genomic position columns must be integers,
    ## NOT numeric (which can cause issues when tabix-indexing).
    bp_cols <- unique(c(start_col, end_col))
    sumstats_dt[, (bp_cols):=lapply(.SD, as.integer), 
                .SDcols=bp_cols] 
    ### setorderv is much more efficient than dplyr::arrange
    data.table::setorderv(x = sumstats_dt, 
                          cols = unique(c(chr_col, start_col, end_col)))
    return(sumstats_dt) 
}