#' Ensure all SNPs have N less than X std dev below mean
#' 
#'Incase some SNPs were genotyped by a specialized genotyping array and 
#'have substantially more samples than others. These will be removed.
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param N_std numeric The number of standard deviations above the mean a SNP's N is needed to be removed. Default is 5.
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom stats sd
check_n_num <- function(sumstats_dt, path, N_std){
  N = NULL
  col_headers <- names(sumstats_dt)
  if("N" %in% col_headers && N_std>0){
    mean_N <- mean(sumstats_dt$N)
    sd_N <- stats::sd(sumstats_dt$N)
    num_bad_ids <- nrow(sumstats_dt[N>((N_std*sd_N)+mean_N),])
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs have N values ",
                    N_std," standard deviations above the mean",
                    " and will be removed")
      message(msg)
      sumstats_dt <- sumstats_dt[N<=((N_std*sd_N)+mean_N),]
    }
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
