granges_style <- function(gr,
                          style = c("NCBI", "UCSC")) {
    suppressWarnings(
        GenomeInfoDb::seqlevelsStyle(gr) <- style[1]
    )
    return(gr)
}
