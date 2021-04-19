#' Remove SNPs with strand-ambiguous alleles
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38). Default is GRCh37.
#' @param strand_ambig_filter Binary Should SNPs with strand-ambiguous alleles be removed. Default is FALSE
#' @return The modified sumstats_file
#' @importFrom data.table setkey
#' @importFrom data.table :=
check_strand_ambiguous <- 
  function(sumstats_file, path, ref_genome, strand_ambig_filter){
  CHR = alleles_as_ambig = SNP = A1 = A2 = NULL
  # If SNP present and user specified to remove
  col_headers <- names(sumstats_file)
  if("SNP" %in% col_headers && !isFALSE(strand_ambig_filter)){
    A_T_ambig <- sumstats_file[A1=="A"&A2=="T" | A1=="T"&A2=="A",]$SNP
    C_G_ambig <- sumstats_file[A1=="C"&A2=="G" | A1=="G"&A2=="C",]$SNP
    num_bad_ids <- length(A_T_ambig)+length(C_G_ambig)
    #check for SNPs not on ref genome
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs are strand-ambiguous alleles including",
                      " ",length(A_T_ambig)," A/T and ",length(C_G_ambig),
                      " C/G ambiguous SNPs. These will be removed")
      message(msg)
      rmv_snps <- c(A_T_ambig, C_G_ambig)
      # join using SNP
      data.table::setkey(sumstats_file,SNP)
      #remove strand ambiguous SNPs
      sumstats_file <- sumstats_file[!rmv_snps,]
    }
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
