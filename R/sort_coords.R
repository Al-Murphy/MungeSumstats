#' Sort sum stats
#' 
#' Sort summary statistics table by genomic coordinates.
#' @param sumstats_dt \link[data.table]{data.table} obj of the
#'  summary statistics file for the GWAS.
#' @param sort_coords Whether to sort by coordinates.
#' @param sort_method Method to sort coordinates by: 
#' \itemize{
#' \item{"data.table" (default)}{Uses \link[data.table]{setorderv},
#'  which is must faster than "GenomicRanges" 
#'  but less robust to variations in some sum stats files.}
#' \item{"GenomicRanges"}{Uses \link[GenomicRanges]{sort.GenomicRanges},
#'  which is more robust to variations in sum stats files
#'   but much slower than the "data.table" method.}
#' }
#' @returns Sorted sumstats_dt
#' 
#' @keywords internal
#' @importFrom data.table :=
sort_coords <- function(sumstats_dt,
                        sort_coordinates = TRUE, 
                        sort_method=c("data.table","GenomicRanges")) {
  ### Add this to avoid confusing BiocCheck
  CHR <- NULL
  if (isTRUE(sort_coordinates)) {
    #### Report ####
    sort_method <- sort_method[1]
    messager("Sorting coordinates with",paste0(shQuote(sort_method),".")) 
    ### Double check that X and Y are uppercase
    sumstats_dt[, CHR := gsub("x|23", "X", CHR)]
    sumstats_dt[, CHR := gsub("y", "Y", CHR)]
    sumstats_dt[, CHR := gsub("mt", "MT", CHR)]
    #### Sort ####
    if(sort_method=="data.table"){
      sumstats_dt <- sort_coords_datatable(sumstats_dt = sumstats_dt)   
    } else if (sort_method=="GenomicRanges"){
      sumstats_dt <- sort_coord_genomicranges(sumstats_dt = sumstats_dt)   
    }
    ### Now set CHR back to character to avoid issues 
    # when merging with other dts
    sumstats_dt[,CHR:=as.character(CHR)]
  }
  return(sumstats_dt)
}
