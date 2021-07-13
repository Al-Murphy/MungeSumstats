#' Check if save path is appropriate
#'
#' Returns corrected \code{save_path}, the file type, and the separator.
#' @inheritParams format_sumstats
#' @keywords internal 
check_save_path <- function(save_path,
                            write_vcf=FALSE){
    suffixes <- supported_suffixes()
    if(is.null(save_path)){ 
        save_path <- paste0(tempfile(),".tsv.gz")
        file_type <- "tempfile" 
        sep <- "\t"
    } else {  
        suffix_match <- sapply(suffixes, function(x){grepl(x, tolower(save_path), ignore.case = TRUE)}) 
        if(sum(suffix_match)>0){
            file_type <- names(suffix_match)[suffix_match][1]
            sep <- if(file_type==".csv") "," else "\t"
        } else {
            if(write_vcf==FALSE) { 
                stop("save_path file format not recognized: ",save_path,
                     "\nMust be one of: \n   ", paste(suffixes, collapse = "\n   "))
            } 
        }
    } 
    if(write_vcf){ 
        save_path <- gsub(paste(suffixes,collapse = "|"),".vcf.gz", save_path)
        sep = "\t"
        file_type = "vcf"
    } 
    
    message("Results will be saved to ==> ",save_path) 
    return(list(
        save_path=save_path,
        file_type=file_type,
        sep=sep
    ))
}
