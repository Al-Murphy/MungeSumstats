#' Ensure all SNPs on specified chromosomes are removed
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param rmv_chr vector or character The chromosomes on which the SNPs should be removed. Use NULL if no filtering necessary. Default is X, Y and mitochondrial. 
#' @return The modified sumstats_file
check_chr <- function(sumstats_file, path, rmv_chr){
  CHR = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_file)
  if("CHR" %in% col_headers && !is.null(rmv_chr)){
    #check for chromosomes to be removed
    if(any(rmv_chr %in% unique(sumstats_file$CHR))){
      num_bad_ids <- nrow(sumstats_file[CHR %in% (rmv_chr),])
      msg <- paste0(num_bad_ids, " SNPs are on chromosomes ",
                    paste(rmv_chr,collapse = ", "),
                    " and will be removed")
      message(msg)
      #remove rows on these chromosomes
      sumstats_file <- sumstats_file[!CHR %in% (rmv_chr),]
    }
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
