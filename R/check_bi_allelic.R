#' Remove non-biallelic SNPs
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38). Default is GRCh37.
#' @param bi_allelic_filter Binary Should non-biallelic SNPs be removed. Default is TRUE
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
#' @importFrom Biostrings IUPAC_CODE_MAP
check_bi_allelic <- 
  function(sumstats_dt, path, ref_genome, bi_allelic_filter, rsids){
    CHR = alleles_as_ambig = SNP = NULL
    # If SNP present and user specified to remove
    col_headers <- names(sumstats_dt)
    if("SNP" %in% col_headers && !isFALSE(bi_allelic_filter)){
      message("Checking for bi-allelic SNPs.")
      if(is.null(rsids)){
        rsids <- load_ref_genome_data(data.table::copy(sumstats_dt$SNP), ref_genome)
      }
      #get chars for SNPs not bi/tri allelic or strand ambig from IUPAC_CODE_MAP
      nonambig_IUPAC_CODE_MAP <- 
        names(Biostrings::IUPAC_CODE_MAP[nchar(Biostrings::IUPAC_CODE_MAP)<3])
      #ensure rsids is up-to-date with filtered sumstats_dt
      rsids <- rsids[unique(sumstats_dt$SNP),,nomatch=NULL]
      data.table::setkey(rsids,SNP)
      num_bad_ids <- nrow(rsids[!alleles_as_ambig %in% nonambig_IUPAC_CODE_MAP])
      #check for SNPs not on ref genome
      if(num_bad_ids>0){
        msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs are non-biallelic.",
                      " These will be removed.")
        message(msg)
        # join using SNP
        data.table::setkey(sumstats_dt,SNP)
        keep_snps <- rsids[alleles_as_ambig %in% nonambig_IUPAC_CODE_MAP]$SNP
        #remove rows missing from the reference genome
        sumstats_dt <- sumstats_dt[keep_snps,]
      }
      return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
    }
    else{
      return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
    }
  }
