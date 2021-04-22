#' Ensure all SNPs on specified chromosomes are removed
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param rmv_chr vector or character The chromosomes on which the SNPs should be removed. Use NULL if no filtering necessary. Default is X, Y and mitochondrial. 
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
check_chr <- function(sumstats_dt, path, rmv_chr){
  CHR = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_dt)
  if("CHR" %in% col_headers && !is.null(rmv_chr)){
    #check for chromosomes to be removed
    if(any(rmv_chr %in% unique(sumstats_dt$CHR))){
      num_bad_ids <- nrow(sumstats_dt[CHR %in% (rmv_chr),])
      msg <- paste0(num_bad_ids, " SNPs are on chromosomes ",
                    paste(rmv_chr,collapse = ", "),
                    " and will be removed")
      message(msg)
      #remove rows on these chromosomes
      sumstats_dt <- sumstats_dt[!CHR %in% (rmv_chr),]
    }
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
