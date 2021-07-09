#' Read in file header
#'
#' @inheritParams format_sumstats
#' @inheritParams readLines 
#' @keywords internal 
read_header <- function(path, n=2){
    message("Reading header.")
    header <- readLines(con = path, n = n)
    # path <- "~/Desktop/ewce/MAGMA_Celltyping/ieu-b-2.vcf.gz"
    #Deal with strange, not recognised characters in header like '\' e.g 'xa6\xc2'
    header[[1]] <- iconv(enc2utf8(header[[1]]),sub="byte") 
    return(header)
}