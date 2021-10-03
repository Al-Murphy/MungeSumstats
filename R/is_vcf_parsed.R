is_vcf_parsed <- function(sumstats_file,
                          verbose=TRUE){
    is_parsed <- "PARSED" %in% colnames(sumstats_file)
    if(is_parsed){
        if(verbose){
            message("sumstats_file previously parsed with vcf2df.")     
        }
    } 
    return(is_parsed)
}