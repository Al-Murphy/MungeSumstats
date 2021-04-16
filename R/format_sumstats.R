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
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setcolorder
#' @importFrom utils read.table
#' @import stringr
#' @export
format_sumstats <- function(path,ref_genome="GRCh37", convert_small_p=TRUE,
                              convert_n_int=TRUE, analysis_trait=NULL){
  #Check input parameters
  validate_parameters(path,ref_genome, convert_small_p,
                        convert_n_int, analysis_trait)
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

  # Check 7: check if CHR and BP are missing but SNP is present
  sumstats_file <- check_no_chr_bp(sumstats_file, path, ref_genome)

  # Check 8: check if CHR and BP are present but SNP is missing
  sumstats_file <- check_no_snp(sumstats_file, path, ref_genome)
  
  # Check 9: check if SNP is present but A1 and/or A2 is missing
  sumstats_file <- check_no_allele(sumstats_file, path, ref_genome)

  # Check 10: check that all the vital columns are present
  check_vital_col(sumstats_file)

  # Check 11: check there is at least one signed sumstats column
  check_signed_col(sumstats_file)

  # Check 12: check first three column headers are SNP, CHR, BP (in that order)
  sumstats_file <- check_col_order(sumstats_file, path)

  #Check 13: Keep only rows which have the number of columns expected
  sumstats_file <- check_miss_data(sumstats_file, path)

  # The formatting process can (rarely) result in duplicated columns,
  # i.e. CHR, if CHR:BP is expanded and one already exists... delete duplicates
  # Check 14: check for duplicated columns
  sumstats_file <- check_dup_col(sumstats_file, path)

  #Check 15: check for small P-values (3e-400 or lower)
  sumstats_file <- check_small_p_val(sumstats_file, path, convert_small_p)

  #Check 16: check is N column not all integers, if so round it up
  sumstats_file <- check_n_int(sumstats_file, path, convert_n_int)

  #Check 17: check all rows have SNPs starting with SNP or rs, drop those don't
  sumstats_file <- check_row_snp(sumstats_file, path)

  #Check 18: check all rows for duplicated SNPs, remove any that are
  sumstats_file <- check_dup_snp(sumstats_file, path)

  # Show how the data now looks
  message("Succesfully finished preparing sumstats file, preview:")
  rows_of_data <- c(sumstats_file[1], sumstats_file[2], sumstats_file[2])
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  headers <- paste(col_headers,collapse = "\t")
  message(headers)
  line1_data <- paste(strsplit(sumstats_file[2], "\t")[[1]],collapse = "\t")
  line2_data <- paste(strsplit(sumstats_file[3], "\t")[[1]],collapse = "\t")
  message(line1_data)
  message(line2_data)

  return(tmp) # Returns address of modified file
}




