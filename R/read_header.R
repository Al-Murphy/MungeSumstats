#' Read in file header
#'
#' @return First \code{n} lines of the VCF header
#'
#' @inheritParams format_sumstats
#' @param n integer. The (maximal) number of lines to read. Negative values
#' indicate that one should read up to the end of input on the connection.
#' @param skip_vcf_metadata logical, should VCF metadata be ignored 
#' @inheritParams format_sumstats
#' 
#' @importFrom VariantAnnotation scanVcfHeader
#' @export
#' @examples
#' path <- system.file("extdata", "eduAttainOkbay.txt", 
#'                     package = "MungeSumstats") 
#' header <- read_header(path = path)                    
read_header <- function(path,
                        n = 2L,
                        skip_vcf_metadata = FALSE,
                        nThread = 1) {
    if(is.null(n)) {
        message("Reading entire file.")
        n <- -2L
    } else {
        message("Reading header.")
    } 
    vcf_suffixes <- supported_suffixes(tabular = FALSE,
                                       tabular_compressed = FALSE)
    if (startsWith(path, "https://gwas.mrcieu.ac.uk") | 
        any(endsWith(path, vcf_suffixes))) {
        #### Read VCFs ####
        if (isTRUE(skip_vcf_metadata)) {  
            if(endsWith(path,".bgz")){
                header <- data.table::fread(text = readLines(path,
                                                              n = n+1000),
                                            nrows = n,
                                            skip = "#CHR",
                                            nThread = nThread)
            } else {
                header <- data.table::fread(input = path,
                                            nrows = n,
                                            skip = "#CHR",
                                            nThread = nThread)
            } 
        } else {
            #### Attempt to read in VCF header as well
            header <- readLines(path, n = 1000)
        }
    } else if (endsWith(path,".bgz")){
        #### Read tabix-indexed tabular ####  
        header <- data.table::fread(text=readLines(con = path,
                                                   n = n+1L)) 
    } else if(any(startsWith(path, c("https:","http:")))){
        #### Read generic remote sumstats ####
        header <- data.table::fread(path,nrows = n) 
    } else {
        #### Read tabular #### 
        header <- readLines(con = path,
                            n = n)
        # Deal with strange, not recognised characters in header
        # like '\' e.g 'xa6\xc2'
        header[[1]] <- iconv(enc2utf8(header[[1]]), sub = "byte")
        if(length(header)<2) header <- rep(header,2)
        header <- data.table::fread(text = header)
    }
    return(header)
}
