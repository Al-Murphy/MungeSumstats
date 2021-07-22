#' Check for N column if not present and user wants, impute N based on user's 
#' sample size. **NOTE** this will be the same value for each SNP which is not
#' necessarily correct and may cause issues down the line.
#' 
#' @return \code{list("sumstats_dt"=sumstats_dt)}
#' 
#' @param sumstats_dt data table obj of the summary statistics file for the 
#' GWAS.
#' @param imputation_ind Binary Should a column be added for each imputation 
#' step to show what SNPs have imputed values for differing fields. This 
#' includes a field denoting SNP allele flipping (flipped). **Note** 
#' these columns will be in the formatted summary statistics returned. Default 
#' is FALSE.
#' @param compute_n Whether to impute N. Default of 0 won't impute, any other 
#' integer will be imputed as the N (sample size) for every SNP in the dataset. 
#' **Note** that imputing the sample size for every SNP is not correct and 
#' 
#' @keywords internal
compute_nsize <- function(sumstats_dt,
                         imputation_ind,
                         compute_n=TRUE){   
    ## Set variables to be used in in place data.table functions to NULL 
    ## to avoid confusing BiocCheck.
    Z = BETA = P = IMPUTATION_n = N = NULL
    
    if(compute_n>0){
        col_headers <- names(sumstats_dt)
        if(!"N" %in% col_headers){
            message("Adding N for all SNPs from sample size")
            sumstats_dt[,N:=compute_n]
            #if user wants information, give SNPs where Z-score calculated
            if(imputation_ind)
                sumstats_dt[,IMPUTATION_n:=TRUE]
        }
    } 
    return(list("sumstats_dt"=sumstats_dt))
}
