#' Check that summary statistics from GWAS are in a homogeneous format
#'
#' @return The address for the modified sumstats file
#'
#' @examples
#' #save the Educational Attainment Okbay sumstat file to a temp directory
#' tmp <- tempfile()
#' writeLines(MungeSumstats::eduAttainOkbay,con = tmp)
#' #pass path to format_sumstats
#' reformatted <- MungeSumstats::format_sumstats(tmp,ref_genome="GRCh37")
#' #returned location has the updated summary statistics file
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38). Default is GRCh37.
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
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setcolorder
#' @importFrom utils read.table
#' @export
format_sumstats <- function(path,ref_genome="GRCh37", convert_small_p=TRUE,
                              convert_n_int=TRUE, analysis_trait=NULL,
                              INFO_filter=0.9, N_std=5, 
                              rmv_chr=c("X","Y","MT"),on_ref_genome=TRUE,
                              strand_ambig_filter=FALSE, allele_flip_check=TRUE,
                              bi_allelic_filter=TRUE
                            ){
  #Avoid reloading ref genome every time, save it to this parent environment
  #after being made once - speed up code
  rsids = NULL
  #Check input parameters
  validate_parameters(path,ref_genome, convert_small_p, convert_n_int, 
                        analysis_trait, INFO_filter, N_std, rmv_chr,
                        on_ref_genome, strand_ambig_filter, allele_flip_check,
                        bi_allelic_filter)
  # This almost surely modifies the file (since most sumstats from different
  # studies are differently formatted), so it makes more sense to just make a
  # temporary file <tmp>, and return the address of the temp
  sumstats_file <- readLines(path)
  #Deal with strange, not recognised characters in header like '\' e.g 'xa6\xc2'
  sumstats_file[[1]] <- iconv(enc2utf8(sumstats_file[[1]]),sub="byte") 
  tmp <- tempfile()
  writeLines(sumstats_file, con=tmp)
  path <- tmp
  

  # Check 1: Check if the file is in VCF format
  sumstats_file <- check_vcf(sumstats_file, path)

  # Check 2: Ensure that tabs separate rows
  sumstats_file <- check_tab_delimited(sumstats_file)

  # Check 3:Standardise headers for all OS
  sumstats_file <-
    standardise_sumstats_column_headers_crossplatform(sumstats_file, path)

  # Check 4: Check if multiple models used or multiple traits tested in GWAS
  sumstats_file <-  check_multi_gwas(sumstats_file, path, analysis_trait)

  col_headers <- sumstats_file[1]
  col_headers <- strsplit(col_headers, "\t")[[1]]

  # Series of checks if CHR or BP columns aren't present
  if(!sum(c("CHR","BP") %in% col_headers)==2){
    msg <- paste0("Summary statistics file does not have obvious CHR/BP colum",
                  "ns. Checking to see if they are joined in another column")
    message(msg)

    #Check 5: check if CHR:BP:A2:A1 merged to 1 column
    sumstats_file <- check_four_step_col(sumstats_file, path) 

    # Check 6: check if there is a column of data with CHR:BP format
    sumstats_file <- check_two_step_col(sumstats_file, path) 

    # Re-standardise in case the joined column headers were unusual
    sumstats_file <-
      standardise_sumstats_column_headers_crossplatform(sumstats_file, path)
  }
  
  #ALL CHECKS BELOW HERE ARE MADE THROUGH DT
  sumstats_file <- data.table::fread(path)

  # Check 7: check if CHR and BP are missing but SNP is present
  sumstats_file <- check_no_chr_bp(sumstats_file, path, ref_genome,
                                    rsids)
  
  # Check 8: check if CHR and BP are present but SNP is missing
  sumstats_file <- check_no_snp(sumstats_file, path, ref_genome)
  
  # Check 9: check if SNP is present but A1 and/or A2 is missing
  sumstats_file <- check_no_allele(sumstats_file, path, ref_genome,
                                    rsids)
  
  # Check 10: check that all the vital columns are present
  check_vital_col(sumstats_file)
  
  # Check 11: check there is at least one signed sumstats column
  check_signed_col(sumstats_file)
  
  # Check 12: check for allele flipping
  sumstats_file <- check_allele_flip(sumstats_file, path, ref_genome,
                                     rsids,allele_flip_check)

  # Check 13: check first three column headers are SNP, CHR, BP (in that order)
  sumstats_file <- check_col_order(sumstats_file, path)

  #Check 14: Keep only rows which have the number of columns expected
  sumstats_file <- check_miss_data(sumstats_file, path)

  # The formatting process can (rarely) result in duplicated columns,
  # i.e. CHR, if CHR:BP is expanded and one already exists... delete duplicates
  # Check 15: check for duplicated columns
  sumstats_file <- check_dup_col(sumstats_file, path)
  
  #Check 16: check for small P-values (3e-400 or lower)
  sumstats_file <- check_small_p_val(sumstats_file, path, convert_small_p)

  #Check 17: check is N column not all integers, if so round it up
  sumstats_file <- check_n_int(sumstats_file, path, convert_n_int)

  #Check 18: check all rows have SNPs starting with SNP or rs, drop those don't
  sumstats_file <- check_row_snp(sumstats_file, path)

  #Check 19: check all rows for duplicated SNPs, remove any that are
  sumstats_file <- check_dup_snp(sumstats_file, path)
  
  #Check 20: check all rows for duplicated BPs, remove any that are
  sumstats_file <- check_dup_bp(sumstats_file, path)
  
  #Check 21: check for low INFO scores
  sumstats_file <- check_info_score(sumstats_file, path, INFO_filter)
  
  #Check 22: check for N > X std dev above mean
  sumstats_file <- check_n_num(sumstats_file, path, N_std)
  
  #Check 23: check that no snps are on specific chromosomes
  sumstats_file <- check_chr(sumstats_file, path, rmv_chr)
  
  #Check 24: check that all snps are present on reference genome
  sumstats_file <- 
    check_on_ref_genome(sumstats_file, path, ref_genome, on_ref_genome,rsids)
  
  #Check 25: check that all snps are not strand ambiguous
  sumstats_file <- 
    check_strand_ambiguous(sumstats_file, path, ref_genome, strand_ambig_filter)
  
  #Check 26: check for non-biallelic SNPS
  sumstats_file <- 
    check_bi_allelic(sumstats_file, path, ref_genome, bi_allelic_filter,
                           rsids)
  
  #WRITE DT DATA TO PATH
  data.table::fwrite(sumstats_file,file=path,sep="\t")
  #read in the first few lines to preview
  sumstats_file <- readLines(path,5L)
  
  # Show how the data now looks
  message("Succesfully finished preparing sumstats file, preview:")
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  headers <- paste(col_headers,collapse = "\t")
  message(headers)
  line1_data <- paste(strsplit(sumstats_file[2], "\t")[[1]],collapse = "\t")
  line2_data <- paste(strsplit(sumstats_file[3], "\t")[[1]],collapse = "\t")
  message(line1_data)
  message(line2_data)

  return(tmp) # Returns address of modified file
}




