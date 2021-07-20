#' Ensure that the input parameters are logical
#' 
#' @inheritParams format_sumstats
#' @return No return
#' @keywords internal
validate_parameters <- function(path,
                                ref_genome, 
                                convert_small_p,
                                compute_z,
                                convert_n_int, 
                                analysis_trait, 
                                INFO_filter,
                                pos_se,
                                effect_columns_nonzero,
                                N_std, 
                                N_dropNA,
                                rmv_chr, 
                                on_ref_genome, 
                                strand_ambig_filter,
                                allele_flip_check,
                                allele_flip_drop,
                                allele_flip_z,
                                bi_allelic_filter,
                                snp_ids_are_rs_ids,
                                write_vcf,
                                return_format,
                                ldsc_format){
  # Checking if the file exists should happen first
  if (!file.exists(path) && !startsWith(path,"https://gwas.mrcieu.ac.uk"))
    stop("Path to GWAS sumstats is not valid")
  
  gen_msg <- paste0("The chosen genome build must be one of GRCh37 or GRCh38 ",
                      "or left as null so the genome build will be inferred ",
                      "from the data.")
  #Check genome build is valid option
  if(!is.null(ref_genome) && !(toupper(ref_genome) %in% c("GRCH37","GRCH38")))
    stop(gen_msg)
  
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
  
  if(is.null(ref_genome)){
    if(!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37", quietly = TRUE))
      stop(GRCH37_msg1)
    if(!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE))
      stop(GRCH37_msg2)
    if(!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38", quietly = TRUE))
      stop(GRCH38_msg1)
    if(!requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE))
      stop(GRCH38_msg2)
  }


  #Check binary values  
  if(!is.logical(N_dropNA))
    stop("N_dropNA must be either TRUE or FALSE")
  if(!is.logical(convert_small_p))
    stop("convert_small_p must be either TRUE or FALSE")
  if(!is.logical(compute_z))
    stop("compute_z must be either TRUE or FALSE") 
  if(!is.logical(convert_n_int))
    stop("convert_n_int must be either TRUE or FALSE")
  if(!is.logical(on_ref_genome))
    stop("on_ref_genome must be either TRUE or FALSE")
  if(!is.logical(strand_ambig_filter))
    stop("strand_ambig_filter must be either TRUE or FALSE")
  if(!is.logical(allele_flip_check))
    stop("allele_flip_check must be either TRUE or FALSE")
  if(!is.logical(allele_flip_drop))
    stop("allele_flip_drop must be either TRUE or FALSE")
  if(!is.logical(allele_flip_z))
    stop("allele_flip_z must be either TRUE or FALSE")
  if(!is.logical(bi_allelic_filter))
    stop("bi_allelic_filter must be either TRUE or FALSE")
  if(!is.logical(snp_ids_are_rs_ids))
    stop("snp_ids_are_rs_ids must be either TRUE or FALSE")
  if(!is.logical(write_vcf))
    stop("write_vcf must be either TRUE or FALSE")
  if(!is.logical(ldsc_format))
    stop("ldsc_format must be either TRUE or FALSE")
  if(!is.logical(pos_se))
    stop("pos_se must be either TRUE or FALSE")
  if(!is.logical(effect_columns_nonzero))
    stop("effect_columns_nonzero must be either TRUE or FALSE")
  
  
  effect_columns_nonzero
  
  
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
  #check return_format
  rf_msg <- paste0("MungeSumstats can only output as data.table, GRanges or ",
                    "VRanges objects, not ",return_format,
                    ". Check return_format")
  if(!(tolower(return_format) %in% c("vr","vranges","gr","granges",
                                      "genomicranges","data.table"))) 
    stop(rf_msg)


}
