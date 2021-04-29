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
#' @importFrom data.table rbindlist
#' @importFrom data.table setcolorder
check_on_ref_genome <- 
  function(sumstats_dt, path, ref_genome, on_ref_genome, rsids){
  CHR = SNP = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_dt)
  if("SNP" %in% col_headers && !isFALSE(on_ref_genome)){
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_dt$SNP), ref_genome)
    }
    #ensure rsids is up-to-date with filtered sumstats_dt
    rsids <- rsids[sumstats_dt$SNP,,nomatch=NULL]
    data.table::setkey(rsids,SNP)
    num_bad_ids <- length(sumstats_dt$SNP) - length(rsids$SNP)
    #check for SNPs not on ref genome
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs are not on the reference genome. ",
                    " These will be corrected from the reference genome.")
      message(msg)
      # join using SNP
      data.table::setkey(sumstats_dt,SNP)
      #if the dataset has CHR & BP, try impute the correct ones
      if(sum(c("CHR","BP") %in% col_headers)==2){
        bad_snp <- sumstats_dt[!rsids$SNP,]
        #remove snp column and pass to function to impute snp
        bad_snp <- bad_snp[,SNP:=NULL]
        corrected_snp <- 
          check_no_snp(bad_snp, path=tempfile(), ref_genome, verbose=FALSE)
        corrected_snp <- corrected_snp$sumstats_dt 
        #make sure columns in correct order
        data.table::setcolorder(corrected_snp,names(sumstats_dt))
        #remove rows missing from the reference genome and combine
        sumstats_dt <- 
          data.table::rbindlist(list(sumstats_dt[rsids$SNP,],corrected_snp))
      }
      else{
        msg <- paste0(num_bad_ids, " SNPs are not on the reference genome. ",
                    " These will be removed")
        message(msg)
        
        #remove rows missing from the reference genome
        sumstats_dt <- sumstats_dt[rsids$SNP,]
      }
    }
  }
  return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
}
