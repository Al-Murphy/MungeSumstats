#' Determine sum stats file type and reda them into memory
#'
#' @inheritParams format_sumstats
#' @keywords internal 
read_sumstats <- function(path,   
                            nThread=1){
    header <- read_header(path = path)
    is_vcf <- check_vcf(header = header)  
    if(is_vcf){
        message("Importing.")
        sumstats_file <- read_vcf(path = path, nThread = nThread)
    } else {
        is_tabular <- check_tabular(header = header) 
        if(is_tabular){
            message("Importing.")
            sumstats_file <- data.table::fread(path, nThread = nThread)
        } else {
            suffixes <- supported_suffixes()
            stop("Unrecognized file format.\n",
                 "Must be one of: \n   ", paste(suffixes, collapse = "\n   "))  
        } 
    } 
    return(sumstats_file) 
}
