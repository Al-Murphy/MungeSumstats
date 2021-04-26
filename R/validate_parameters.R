#' Ensure that the input parameters are logical
#'
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param convert_small_p Binary, should p-values < 5e-324 be converted to 0? Small p-values pass the R limit and can cause errors with LDSC/MAGMA and should be converted. Default is TRUE.
#' @param convert_n_int Binary, if N (the number of samples) is not an integer, should this be rounded? Default is TRUE.
#' @param analysis_trait If multiple traits were studied, name of the trait for analysis from the GWAS. Default is NULL
#' @param INFO_filter numeric The minimum value permissible of the imputation information score (if present in sumstatsfile). Default 0.9
#' @param N_std numeric The number of standard deviations above the mean a SNP's N is needed to be removed. Default is 5.
#' @param rmv_chr vector or character The chromosomes on which the SNPs should be removed. Use NULL if no filtering necessary. Default is X, Y and mitochondrial. 
#' @param on_ref_genome Binary Should a check take place that all SNPs are on the reference genome by SNP ID. Default is TRUE 
#' @param strand_ambig_filter Binary Should SNPs with strand-ambiguous alleles be removed. Default is FALSE
#' @param allele_flip_check Binary Should the allele columns be chacked against reference genome to infer if flipping is necessary. Default is TRUE
#' @param bi_allelic_filter Binary Should non-biallelic SNPs be removed. Default is TRUE
#' @return No return
#' @keywords internal
validate_parameters <- function(path,ref_genome, convert_small_p,
                            convert_n_int, analysis_trait, INFO_filter,
                            N_std, rmv_chr, on_ref_genome, strand_ambig_filter,
                            allele_flip_check, bi_allelic_filter){
  # Checking if the file exists should happen first
  if (!file.exists(path))
    stop("Path to GWAS sumstats is not valid")

  #Check genome build is valid option
  if(!(toupper(ref_genome) %in% c("GRCH37","GRCH38")))
    stop("The chosen genome build must be one of GRCh37 or GRCh38")
  
  #checks for installed packages
  GRCH37_msg1 <- paste0("Install 'SNPlocs.Hsapiens.dbSNP144.GRCh37' to use ",
                          "'GRCh37' as 'ref_genome'")
  GRCH37_msg2 <- paste0("Install 'BSgenome.Hsapiens.1000genomes.hs37d5' to ",
                          "use 'GRCh37' as 'ref_genome'")
  if(toupper(ref_genome)=="GRCH37" && 
      !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37", quietly = TRUE))
    stop(GRCH37_msg1)
  if(toupper(ref_genome)=="GRCH37" && 
      !requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE))
    stop(GRCH37_msg2)
  
  GRCH38_msg1 <- paste0("Install 'SNPlocs.Hsapiens.dbSNP144.GRCh38' to use ",
                          "'GRCh38' as 'ref_genome'")
  GRCH38_msg2 <- paste0("Install 'BSgenome.Hsapiens.NCBI.GRCh38' to ",
                          "use 'GRCh38' as 'ref_genome'")
  if(toupper(ref_genome)=="GRCH38" && 
     !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38", quietly = TRUE))
    stop(GRCH38_msg1)
  if(toupper(ref_genome)=="GRCH38" && 
     !requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE))
    stop(GRCH38_msg2)


  #Check binary values
  if(!is.logical(convert_small_p))
    stop("convert_small_p must be either TRUE or FALSE")
  if(!is.logical(convert_n_int))
    stop("convert_n_int must be either TRUE or FALSE")
  if(!is.logical(on_ref_genome))
    stop("on_ref_genome must be either TRUE or FALSE")
  if(!is.logical(strand_ambig_filter))
    stop("strand_ambig_filter must be either TRUE or FALSE")
  if(!is.logical(allele_flip_check))
    stop("allele_flip_check must be either TRUE or FALSE")
  if(!is.logical(bi_allelic_filter))
    stop("bi_allelic_filter must be either TRUE or FALSE")
  
  
  #Check numeric
  if(!is.numeric(INFO_filter))
    stop("INFO_filter must be numeric")
  if(!is.numeric(N_std))
    stop("N_std must be numeric")
  
  #Check rmv_chr choices all valid chromosomes
  chrs <- c(as.character(seq_len(22)),"X","Y","MT")
  chr_msg <- 
    paste0("rmv_chr choices must be one/or more of: \n",
            paste(chrs, collapse = ", "))
  if(!is.null(rmv_chr))
    if(!all(rmv_chr %in% chrs))
      stop(chr_msg)

}
