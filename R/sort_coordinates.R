#' Sort sum stats by genomic coordinates
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param sort_coords Whether to sort by coordinates.
#' @keywords internal 
#' @importFrom dplyr %>% arrange 
sort_coords <- function(sumstats_dt,
                        sort_coordinates=TRUE){
    if(sort_coordinates){
        message("Sorting coordinates")
        sumstats_sorted <- sumstats_dt %>%
            dplyr::arrange(CHR, POS)
        return(sumstats_sorted)
    } else { return(sumstats_dt) }  
}