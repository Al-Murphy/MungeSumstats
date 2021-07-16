#' Check if the inputted file is in VCF format
#'
#' @param header Header of the GWAS summary statistics file. 
#' @return Whether the file is vcf or not
#' @keywords internal 
#' @importFrom methods is
check_vcf <- function(header){
  if(is(header,"VCFHeader")) return(TRUE)
  P = LP = INFO = NULL
  #if the file is a VCF, first line will look like: ##fileformat=VCFv4.2
  first_line <- header[[1]]
  file_type <- gsub("^##fileformat=","",first_line)
  is_vcf <- length(grep("^vcf",tolower(file_type)))==1
  if(is_vcf){ 
    message("VCF format detected.",
            "This will be converted to a standardised table format.")
  }
  return(is_vcf)
}
