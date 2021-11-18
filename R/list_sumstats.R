#' List munged summary statistics
#' 
#' Searches for and lists local GWAS summary statistics files munged by 
#' \link[MungeSumstats]{format_sumstats} or 
#' \link[MungeSumstats]{import_sumstats}.
#' 
#' @param save_dir Top-level directory to recursively search 
#' for summary statistics files within.
#' @param pattern Regex pattern to search for files with.
#' @param ids_from_file Try to extract dataset IDs from file names.
#' If \code{FALSE}, will infer IDs from the directory names instead.
#' @param verbose Print messages.
#' 
#' @returns Named vector of summary stats paths.
#' 
#' @examples
#' save_dir <- system.file("extdata",package = "MungeSumstats")
#' munged_files <- MungeSumstats::list_sumstats(save_dir = save_dir)
#' 
#' @export
#' @importFrom stats setNames
list_sumstats <- function(save_dir = getwd(),
                          pattern = "*.tsv.gz$",
                          ids_from_file = TRUE,
                          verbose = TRUE){ 
    munged_files <- list.files(path = save_dir,
                               pattern = pattern,
                               full.names = TRUE, 
                               recursive = TRUE)
    #### Exclude log files ####
    munged_files <- munged_files[basename(dirname(munged_files))!="logs"] 
    #### infer IDs ####
    if(ids_from_file){
        ids <- gsub(paste(supported_suffixes(),collapse = "|"), "",
                    basename(munged_files))
    }else {
        ids <- basename(dirname(munged_files))
    }
    #### Report ####
    messager(length(munged_files),"file(s) found.",v=verbose)
    return(stats::setNames(munged_files,ids))
}

