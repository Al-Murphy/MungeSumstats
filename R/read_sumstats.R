#' Determine summary statistics file type and read them into memory
#'
#' @return \code{data.table} of formatted summary statistics
#'
#' @param nrows integer. The (maximal) number of lines to read.
#' If \code{Inf}, will read in all rows.
#' @param standardise_headers Standardise headers first.
#' @inheritParams format_sumstats
#' @inheritParams standardise_header
#' @inheritParams vcf2df
#' @inheritParams read_vcf
#' 
#' @export
#' @importFrom data.table fread as.data.table setkeyv
#' @examples
#' path <- system.file("extdata", "eduAttainOkbay.txt",
#'     package = "MungeSumstats"
#' )
#' eduAttainOkbay <- read_sumstats(path = path)
read_sumstats <- function(path, 
                          nrows = Inf,
                          standardise_headers = FALSE,
                          samples = 1,
                          add_sample_names = FALSE,
                          nThread = 1,
                          mapping_file = sumstatsColHeaders) {
    
    if (is.data.frame(path)) {
        message("Summary statistics passed as R object.")
        sumstats_file <- data.table::as.data.table(path)
        if (!is.infinite(nrows)) {
            sumstats_file <- sumstats_file[seq(1, nrows), ]
        }
    } else {
        vcf_suffixes <- supported_suffixes(tabular = FALSE, 
                                           tabular_compressed = FALSE)
        is_vcf <- grepl(paste(vcf_suffixes,collapse = "|"), path) 
        if (isTRUE(is_vcf)) {
            sumstats_file <- read_vcf(path = path,
                                      use_params = TRUE,
                                      samples = samples,
                                      nThread = nThread)
        } else {
            #### Check if tabular 1: infer from file name ####
            tab_suffixes <- supported_suffixes(vcf = FALSE, 
                                               vcf_compressed = FALSE)
            is_tabular <- grepl(paste(tab_suffixes,collapse = "|"), path) 
            #### Check if tabular 2: infer from data ####
            if(isFALSE(is_tabular)){
                header <- read_header(path = path)
                is_tabular <- check_tabular(header = header)
            }
            #### Process tabular ####
            if (isTRUE(is_tabular)) {
                if(endsWith(path,".bgz")){
                    message("Importing tabular bgz file: ", path)
                    sumstats_file <- data.table::fread(
                        text = readLines(con = path),
                        nThread = nThread,
                        nrows = nrows)
                }else {
                    message("Importing tabular file: ", path)
                    sumstats_file <- data.table::fread(
                        path,
                        nThread = nThread,
                        nrows = nrows)
                }
               
            } else {
                suffixes <- supported_suffixes()
                stop(
                    "Unrecognized file format.\n",
                    "Must be one of: \n   ",
                    paste(suffixes, collapse = "\n   ")
                )
            }
        }
    }
    #### Drop empty cols ####
    sumstats_file <- remove_empty_cols(sumstats_dt = sumstats_file) 
    #### Standardise colnames ####
    if (isTRUE(standardise_headers)) {
        CHR <- NULL;
        sumstats_file <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_file,
                mapping_file = mapping_file
            )[["sumstats_dt"]]
        #### Ensure CHR is character ####
        if("CHR" %in% names(sumstats_file)) {
            sumstats_file[,CHR:=as.character(CHR)]
        }
        #### Ensure SNP is the key ####
        if("SNP" %in% names(sumstats_file)) {
            data.table::setkeyv(sumstats_file, cols = "SNP")
        } 
    }
    return(sumstats_file)
}
