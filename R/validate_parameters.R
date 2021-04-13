#' Ensure that the input parameters are logical
#'
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param convert_small_p Binary, should p-values < 5e-324 be converted to 0? Small p-values pass the R limit and can cause errors with LDSC/MAGMA and should be converted. Default is TRUE.
#' @param convert_n_int Binary, if N (the number of samples) is not an integer, should this be rounded? Default is TRUE.
#' @param analysis_trait If multiple traits were studied, name of the trait for analysis from the GWAS. Default is NULL
#' @return No return
#'

validate_parameters <- function(path,ref_genome, convert_small_p,
                            convert_n_int, analysis_trait){
  # Checking if the file exists should happen first
  if (!file.exists(path))
    stop("Path to GWAS sumstats is not valid")

  #Check genome build is valid option
  if(!(toupper(ref_genome) %in% c("GRCH37","GRCH38")))
    stop("The chosen genome build must be one of GRCh37 or GRCh38")

  #Check binary values
  if(!is.logical(convert_small_p))
    stop("convert_small_p must be either TRUE or FALSE")
  if(!is.logical(convert_n_int))
    stop("convert_n_int must be either TRUE or FALSE")

}
