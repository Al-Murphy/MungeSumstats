#' Sort sum stats by genomic coordinates
#'
#' @return Sorted sumstats_dt
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param sort_coords Whether to sort by coordinates.
#' @param make_ordered Make CHR into an ordered factor to ensure 
#' they go from 1-22, X, Y.
#' @keywords internal
#' @importFrom data.table setorderv
sort_coords <- function(sumstats_dt,
                        sort_coordinates = TRUE) {
    ### Add this to avoid confusing BiocCheck
    CHR <- NULL

    if (sort_coordinates) {
        message("Sorting coordinates.")
        chr_order <- c(seq_len(22), "X", "Y")
        ### Double check that X and Y are uppercase
        sumstats_dt[, CHR := gsub("x|23", "X", CHR)]
        sumstats_dt[, CHR := gsub("y", "Y", CHR)]
        ### Turn CHR into an ordered factor to account for X and Y chroms
        sumstats_dt[, CHR := factor(CHR, levels = chr_order, ordered = TRUE)]
        ### setorderv is much more efficient than dplyr::arrange
        data.table::setorderv(sumstats_dt, c("CHR", "BP"))
        ### Now set CHR back to character to avoid issues 
        # when merging with other dts
        sumstats_dt[, CHR := as.character(CHR)]
    }
    return(sumstats_dt)
}
