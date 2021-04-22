#' Ensure all SNPs are on the reference genome
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38). Default is GRCh37.
#' @param on_ref_genome Binary Should a check take place that all SNPs are on the reference genome by SNP ID. Default is TRUE
#' @param rsids datatable of snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics data table object
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' }
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table copy
check_on_ref_genome <- 
  function(sumstats_dt, path, ref_genome, on_ref_genome, rsids){
  CHR = SNP = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_dt)
  if("SNP" %in% col_headers && !isFALSE(on_ref_genome)){
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_dt$SNP), ref_genome)
    }
    num_bad_ids <- length(sumstats_dt$SNP) - length(rsids$SNP)
    #check for SNPs not on ref genome
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs are not on the reference genome. ",
                    " These will be removed")
      message(msg)
      # join using SNP
      data.table::setkey(sumstats_dt,SNP)
      #remove rows missing from the reference genome
      sumstats_dt <- sumstats_dt[rsids$SNP,]
    }
  }
  return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
}
