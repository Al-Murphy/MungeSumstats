#' Check that summary statistics from GWAS are in a homogeneous format
#'
#' @return The address for the modified sumstats file
#'
#' @examples
#' #format_sumstats(path)
#' @param path Filepath for the summary statistics file to be formatted
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setcolorder
#' @importFrom utils read.table
#' @import stringr
#' @export
format_sumstats <- function(path){
  #Standardise the format of summary statistics from GWAS
  #Homogenise the format of summary statistics from GWAS
  
  #All issues with format_sumstats_for_magma here
  
  
  #use A1/A2 to confirm that it's the same SNP being discussed.... 
  #e.g. if it's CHR 2 bp 234234.... the SNP can be a->c or a->g.... 
  #at the moment I'm not checking that, I just find any SNP at that location
  
  #drop the SNP column from the one in the vignette though to check it works as expected
  
  #also this part of the function.... currently it's commented out... would be great if 
  #you could implement it properly and check it works
  
  #summary statistics from vignette, just remove SNP column so it has to be inferred
  
  
  #test cases:
  #1. column of data with CHR:BP:A2:A1 format - fourStepCol
  #20016_irnt.gwas.imputed_v3.both_sexes.tsv
  
  #Needed test cases:
  # Ensure that tabs separate rows
  
  # Checking if the file exists should happen first
  if (!file.exists(path))
    stop("Path to GWAS sumstats is not valid")
  
  # This almost surely modifies the file (since most sumstats from different 
  # studies are differently formatted), so it makes more sense to just make a 
  # temporary file <tmp>, and return the address of the temp
  sumstats_file <- readLines(path)
  tmp <- tempfile()
  writeLines(sumstats_file, con=tmp)
  path <- tmp
  
  # Check 1: Ensure that tabs separate rows
  sumstats_file <- check_tab_delimited(sumstats_file) 
  
  #Standardise headers for all OS
  sumstats_file[1] <- 
    standardise_sumstats_column_headers_crossplatform(sumstats_file[1])
  col_headers <- sumstats_file[1]
  col_headers <- strsplit(col_headers, "\t")[[1]]
  
  # Series of checks if CHR or BP columns aren't present
  if(!sum(c("SNP","BP") %in% col_headers)==2){
    # - UKBB data from Ben Neale has a Variant column with CHR:POS:REF:ALT where 
    # ALT allele is the effect allele in the model [NB: the ALT allele is NOT always the minor allele]
    # -- For input to LDSC, A1 is effect allele, A2 is non-effect allele
    # - DIAGRAM diabetes data has a Chr:Position column with CHR:BP
    # - BMI adjusted for smoking has markername with CHR:BP (with the chromosome name having 'chr' preceeding)
    # - Agression [EAGLE] just doesn't have any CHR or BP data
    print(paste0("Summary statistics file does not have obvious CHR or BP colu",
                  "mns. Checking to see if they are joined in another column"))
    
    #Check 2: check if CHR:BP:A2:A1 merged to 1 column
    sumstats_file <- check_four_step_col(sumstats_file) 
    
    # Check 3: check if there is a column of data with CHR:BP format
    sumstats_file <- check_two_step_col(sumstats_file) 
    
    # Restandardise in case the joined column headers were unusual
    sumstats_file[1] <- 
      standardise_sumstats_column_headers_crossplatform(sumstats_file[1])
  }
  
  #TODO DON"T KNOW IF THESE LINES ARE NECESSARY??---------
  # If SNP is present... BUT not CHR or BP then need to find the relevant locations
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]) 
  col_headers = strsplit(rows_of_data[1], "\t")[[1]] 
  writeLines(sumstats_file, con = path)
  
  # Check 4: check if CHR and BP are missing but SNP is present
  sumstats_file <- check_no_chr_bp(sumstats_file)
  
  # Check 5: check if CHR and BP are present but SNP is missing
  sumstats_file <- check_no_snp(sumstats_file)
  
  # Check 6: check that all the vital columns are present
  check_vital_col(sumstats_file)
  
  # Check 7: check there is at least one signed sumstats column
  check_signed_col(sumstats_file)
  
  # Check 8: check first three column headers are SNP, CHR, BP (in that order)
  sumstats_file <- check_col_order(sumstats_file, path)
  
  # The formatting process can (rarely) result in duplicated columns, 
  # i.e. CHR, if CHR:BP is expanded and one already exists... delete duplicates
  # Check 9: check for duplicated columns
  sumstats_file <- check_dup_col(sumstats_file, path)
  
  #Check 10: check for small P-values (3e-400 or lower)
  sumstats_file <- check_small_p_val(sumstats_file, path)
  
  #Check 11: check is N column not all integers, if so round it up
  sumstats_file <- check_n_int(sumstats_file, path)
  
  #Check 12: check all rows have SNPs starting with SNP or rs, drop those don't 
  sumstats_file <- check_row_snp(sumstats_file, path)
  
  #TODO - Ask Nathan, I don't think this is necessary - I think DT will throw an error
  # Keep only rows which have the number of columns expected (I've commmented this out because it takes forever to run... but presumably I wrote it for a reason!)
  # print("Keeping only rows which have the number of columns expected.")
  # sumstats_file <- readLines(path); expected_number_of_columns <- length(strsplit(sumstats_file[1],"\t")[[1]]); good_ones <- sumstats_file[1]
  # for (line in sumstats_file) {
  #   if (line == sumstats_file[1]) {next} # Skip header
  #   if ( length(strsplit(line,"\t")[[1]]) == expected_number_of_columns ) {good_ones <- c(good_ones, line)} # This adds every line that has expected number of columns to a temporary list
  # }
  # writeLines(text=good_ones, con = path)
  
  
  #Check 13: check all rows for duplicated SNPs, remove any that are 
  sumstats_file <- check_dup_snp(sumstats_file, path)
  
  # Show how the data now looks
  print("Succesfully finished preparing sumstats file:")
  print("Header of file:")
  con <- file(path,"r")
  rows_of_data <- readLines(con,n=2)
  close(con)
  print(rows_of_data)
  
  return(tmp) # Returns address of modified file
}



#TODO
#check.small.p.val
#check.no.chr.bp
#check.no.snp
#check.n.int
#check.row.snp
#check.dup.snp
#Ask Nathan



