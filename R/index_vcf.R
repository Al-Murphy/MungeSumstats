#' Tabix-index file: VCF
#' 
#' Convert summary stats file to tabix format
#'
#' @source Borrowed function from 
#' \href{https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R}{
#' echotabix}.
#' 
#' @param path Path to VCF.
#' @param verbose Print messages.
#' 
#' @return Path to tabix-indexed tabular file
#'
#' @family tabix
#' @keywords internal
#' @importFrom methods is
#' @importFrom R.utils gunzip
#' @importFrom VariantAnnotation indexVcf
#' @importFrom Rsamtools bgzip
#' @examples 
#' eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
#'                                  package = "MungeSumstats")
#' sumstats_dt <- data.table::fread(eduAttainOkbayPth, nThread = 1)
#' sumstats_dt <- 
#' MungeSumstats:::standardise_sumstats_column_headers_crossplatform(
#'     sumstats_dt = sumstats_dt)$sumstats_dt
#' sumstats_dt <- MungeSumstats:::sort_coords(sumstats_dt = sumstats_dt)
#' path <- tempfile(fileext = ".tsv")
#' MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt, save_path = path)
#'     
#' indexed_file <- MungeSumstats::index_tabular(path = path)
index_vcf <- function(path,
                      verbose=TRUE){ 
    #### When remote ####
    if(!file.exists(path) && 
       any(startsWith(tolower(path),c("http","ftp")))){  
        messager("File inferred to be remote.",v=verbose)
        return(path)
        #### When already tabix ####
    } else if(is_tabix(path = path)){
        messager("File already tabix-indexed.",v=verbose)
        return(path)
        #### When local and non-tabix ####
    } else {
        #### Round 1: assume .gz is bgz-compressed ####
        sfx <- supported_suffixes(tabular = FALSE, 
                                  tabular_compressed = FALSE, 
                                  vcf = FALSE)
        if(!any(endsWith(path,sfx))){
            messager("bgzip-compressing VCF file.",v=verbose)
            path <- Rsamtools::bgzip(file = path, 
                                     dest = sprintf("%s.bgz",
                                                    sub("\\.gz$|\\.bgz$", "",
                                                        path)),
                                     overwrite = TRUE)
        } 
        out <- tryCatch({
            VariantAnnotation::indexVcf(x = path)$path
        }, error = function(e){e}) 
        #### Round 2: assume .gz is NOT bgz-compressed ####
        if(methods::is(out,"error") &&
           grepl("file does not appear to be bgzip",out$message, 
                 ignore.case = TRUE)
           ){
            if(endsWith(tolower(path),".gz")){
                path <- R.utils::gunzip(path, overwrite=TRUE)
            }
            path <- Rsamtools::bgzip(file = path, 
                                     dest = sprintf("%s.bgz",
                                                    sub("\\.gz$|\\.bgz$", "",
                                                        path)),
                                     overwrite = TRUE)
            path <- VariantAnnotation::indexVcf(x = path)$path
            
        }
        return(path)
    } 
}
