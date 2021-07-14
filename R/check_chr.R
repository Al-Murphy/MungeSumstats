#' Ensure all SNPs on specified chromosomes are removed
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param rmv_chr vector or character The chromosomes on which the SNPs should be removed. 
#' Use NULL if no filtering necessary. Default is X, Y and mitochondrial. 
#' @param make_uppercase Make X and Y chromosomes uppercase. 
#' @param rmv_chrPrefix Remove "chr" or "CHR" from chromosome names.
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
check_chr <- function(sumstats_dt, 
                      path, 
                      rmv_chr,
                      make_uppercase=TRUE,
                      rmv_chrPrefix=TRUE){
  CHR = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_dt)
  if("CHR" %in% col_headers && !is.null(rmv_chr)){
    
    ### Sometimes X is labeled as 23
    sumstats_dt[,CHR:=gsub("23","X",CHR)]
    
    #### Remove chr prefix uppercase ####
    if(rmv_chrPrefix){
      message("Removing 'chr' prefix from CHR.")
      sumstats_dt[,CHR:=gsub("chr","",CHR,ignore.case = TRUE)]
      rmv_chr <- gsub("chr","",rmv_chr,ignore.case = TRUE)
    } 
    #### Make all CHR uppercase ####
    if(make_uppercase){
      message("Making X/Y CHR uppercase.")
      sumstats_dt[,CHR:=gsub("x|23","X",CHR)]
      sumstats_dt[,CHR:=gsub("y","Y",CHR)]
    } 
    
    #check for chromosomes to be removed
    ### Standardise chromosomes specified
    rmv_chr <- toupper(rmv_chr)
    if(any(rmv_chr %in% unique(sumstats_dt$CHR))){
      num_bad_ids <- nrow(sumstats_dt[CHR %in% rmv_chr,])
      msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs are on chromosomes ",
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
