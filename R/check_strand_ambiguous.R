#' Remove SNPs with strand-ambiguous alleles
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38). Default is GRCh37.
#' @param strand_ambig_filter Binary Should SNPs with strand-ambiguous alleles be removed. Default is FALSE
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table setkey
#' @importFrom data.table :=
check_strand_ambiguous <- 
  function(sumstats_dt, path, ref_genome, strand_ambig_filter){
  CHR = alleles_as_ambig = SNP = A1 = A2 = NULL
  # If SNP present and user specified to remove
  col_headers <- names(sumstats_dt)
  if("SNP" %in% col_headers && !isFALSE(strand_ambig_filter)){
    A_T_ambig <- sumstats_dt[A1=="A"&A2=="T" | A1=="T"&A2=="A",]$SNP
    C_G_ambig <- sumstats_dt[A1=="C"&A2=="G" | A1=="G"&A2=="C",]$SNP
    num_bad_ids <- length(A_T_ambig)+length(C_G_ambig)
    #check for SNPs not on ref genome
    if(num_bad_ids>0){
      msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs are strand-ambiguous alleles including",
                      " ", formatC(length(A_T_ambig),big.mark = ","),
                    " A/T and ",formatC(length(C_G_ambig),big.mark = ","),
                      " C/G ambiguous SNPs. These will be removed")
      message(msg)
      rmv_snps <- c(A_T_ambig, C_G_ambig)
      # join using SNP
      data.table::setkey(sumstats_dt,SNP)
      #remove strand ambiguous SNPs
      sumstats_dt <- sumstats_dt[!rmv_snps,]
    }
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
