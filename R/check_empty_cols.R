#' Check for empty columns
#'
#' Empty columns contain only ".", NA, or 0
#' @inheritParams format_sumstats 
#' @param sampled_rows First N rows to sample. 
#' Set \code{NULL} to use full \code{sumstats_file}.
#' when determining whether cols are empty. 
#' @return empty_cols
#' @keywords internal 
check_empty_cols <- function(sumstats_file, 
                             sampled_rows=1000){
    if(is.null(sampled_rows)) sampled_rows <- nrow(sumstats_file)
    empty_cols <- sapply(colnames(sumstats_file), function(x){
        (sum(head(sumstats_file, sampled_rows)[[x]]!=".")==0) |
        (sum(!is.na(head(sumstats_file, sampled_rows)[[x]]))==0) |
        (sum(head(sumstats_file, sampled_rows)[[x]]==0)!=0)
    })
    empty_cols <- empty_cols[empty_cols]
    # message(length(empty_cols)," empty columns detected.")
    return(empty_cols)
}