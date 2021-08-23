#' Determine summary statistics file type and read them into memory
#'
#' @return \code{data.table} of formatted summary statistics 
#'
#' @inheritParams format_sumstats
#' @export 
#' @importFrom data.table fread
#' @examples
#'  eduAttainOkbay <- MungeSumstats::read_sumstats(path = system.file("extdata","eduAttainOkbay.txt",
#'                                                                    package="MungeSumstats"))
read_sumstats <- function(path,   
                          nThread=1){
    if(is.data.frame(path)){
        message("Summary statistics passed as R object.")
        sumstats_file<- data.table::as.data.table(path)
    }
    else{
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
                        "Must be one of: \n   ", 
                        paste(suffixes, collapse = "\n   "))  
            }
        }        
    } 
    return(sumstats_file) 
}
