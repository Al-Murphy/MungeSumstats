#' Check for Z-score column
#' 
#' The following ensures that a Z-score column is present.
#' The Z-score formula we used here is a R implementation of the formula
#' used in \href{https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L363}{LDSC's munge_sumstats.py}: 
#' 
#' \code{np.sqrt(chi2.isf(P, 1))} 
#' 
#' The R implementation is adapted from the \code{GenomicSEM::munge} function, 
#' after optimizing for speed using \code{data.table}: 
#' 
#' \code{sumstats_dt[,Z:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))]} 
#' 
#' \emph{NOTE}: \code{compute_z} is set to \code{TRUE} by default to ensure standardisation
#' of the "Z" column (which can be computed differently in different datasets).
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
#' @param compute_z Whether to compute Z-score column from P. Default is FALSE.
#' **Note** that imputing the Z-score for every SNP will not correct be 
#' perfectly correct and may result in a loss of power. This should only be done
#' as a last resort.
#' @param force_new_z When a "Z" column already exists, it will be used by 
#' default. To override and compute a new Z-score column from P set 
#' \code{force_new_z=TRUE}.  
#' @param standardise_headers Run 
#' \code{standardise_sumstats_column_headers_crossplatform} first.  
#' 
#' @keywords internal
#' @importFrom stats qchisq
check_zscore <- function(sumstats_dt,
                         imputation_ind,
                         compute_z=TRUE,
                         force_new_z=FALSE,
                         standardise_headers=FALSE){   
    ## Set variables to be used in in place data.table functions to NULL 
    ## to avoid confusing BiocCheck.
    Z = BETA = P = IMPUTATION_z_score = NULL
    
    if(standardise_headers){
        sumstats_dt <- standardise_sumstats_column_headers_crossplatform(sumstats_dt = sumstats_dt)[["sumstats_dt"]]
    }
    
    if(compute_z){
        # message("Checking Z-score.") 
        col_headers <- names(sumstats_dt)
        if("Z" %in% col_headers && (!force_new_z)){
            message("Keeping existing Z-score column.") 
        } else{
            message("Computing Z-score from P using formula: `sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE)`")
            sumstats_dt[,Z:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))]
            #if user wants information, give SNPs where Z-score calculated
            if(imputation_ind)
                sumstats_dt[,IMPUTATION_z_score:=TRUE]
        }
    } 
    return(list("sumstats_dt"=sumstats_dt))
}
