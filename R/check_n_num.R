#' Ensure all SNPs have N less than X std dev below mean
#' 
#'Incase some SNPs were genotyped by a specialized genotyping array and 
#'have substantially more samples than others. These will be removed.
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param N_std numeric The number of standard deviations above the mean a SNP's N is needed to be removed. Default is 5.
#' @return The modified sumstats_file
#' @importFrom stats sd
check_n_num <- function(sumstats_file, path, N_std){
  N = NULL
  col_headers <- names(sumstats_file)
  if("N" %in% col_headers && N_std>0){
    mean_N <- mean(sumstats_file$N)
    sd_N <- stats::sd(sumstats_file$N)
    num_bad_ids <- nrow(sumstats_file[N>((N_std*sd_N)+mean_N),])
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs have N values ",
                    N_std," standard deviations above the mean",
                    " and will be removed")
      message(msg)
      sumstats_file <- sumstats_file[N<=((N_std*sd_N)+mean_N),]
    }
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
