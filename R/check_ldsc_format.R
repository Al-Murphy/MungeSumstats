#' Ensures that parameters are compatible with LDSC format 
#' 
#' Format summary statistics for direct input to 
#' Linkage Disequilibrium SCore (LDSC) regression without the need
#' to use their \code{munge_sumstats.py} script first. 
#' 
#' \href{https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format}{LDSC documentation}.
#' 
#' @inheritParams format_sumstats  
#' @inheritParams check_zscore 
#' 
#' @return Formatted summary statistics  
#' @source \href{https://github.com/bulik/ldsc}{LDSC GitHub}  
check_ldsc_format <- function(ldsc_format,
                              convert_n_int,
                              allele_flip_check,
                              compute_z
                              ){
    if(ldsc_format){
        message("Ensuring parameters comply with LDSC format.")
        if(!convert_n_int){
            message("Setting `convert_n_int=TRUE` to comply with LDSC format.")
            convert_n_int <- TRUE
        }
        if(!allele_flip_check){
            message("Setting `allele_flip_check=TRUE` to comply with LDSC format.")
            allele_flip_check <- TRUE
        }
        if(!compute_z){
            message("Setting `compute_z=TRUE` to comply with LDSC format.")
            compute_z <- TRUE
        }
    }
    return(list(ldsc_format=ldsc_format,
                convert_n_int=convert_n_int,
                allele_flip_check=allele_flip_check,
                compute_z=compute_z))
}
