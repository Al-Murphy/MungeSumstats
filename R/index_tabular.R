#' Tabix-index file: table
#' 
#' Convert summary stats file to tabix format.
#'
#' @source Borrowed function from 
#' \href{https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R}{
#' echotabix}.
#' 
#' @param path Path to GWAS summary statistics file.
#' @param remove_tmp Remove the temporary uncompressed version of the file 
#' (\emph{.tsv}). 
#' @param verbose Print messages.
#' @inheritParams liftover
#' @inheritParams Rsamtools::bgzip
#' 
#' @return Path to tabix-indexed tabular file
#'
#' @family tabix
#' @export
#' @importFrom R.utils gunzip isGzipped
#' @examples
#' sumstats_dt <- MungeSumstats::formatted_example() 
#' path <- tempfile(fileext = ".tsv")
#' MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt, save_path = path)
#' indexed_file <- MungeSumstats::index_tabular(path = path)
index_tabular <- function(path,
                          chrom_col = "CHR",
                          start_col = "BP",
                          end_col = start_col,
                          overwrite = TRUE,
                          remove_tmp = TRUE,
                          verbose = TRUE) {
    msg <- paste("Converting full summary stats file to",
                 "tabix format for fast querying...")
    messager(msg,v=verbose) 
    #### Make sure input file isn't empty ####
    if(!file.exists(path)){
        stp <- paste("File does not exist: ",path)
        stop(stp)
    }
    if (file.size(path) == 0) {
        msg2 <- paste("Removing empty file:", path)
        messager(msg2,v=verbose)
        out <- file.remove(path)
    }
    #### Ensure it's not already compressed ####
    if(R.utils::isGzipped(path)){
        path <- R.utils::gunzip(path, overwrite = TRUE)
    }
    #### Read header and make dictionary ####
    cdict <- column_dictionary(file_path = path)
    #### Check column exist ####
    if(!chrom_col %in% names(cdict)) stop("chrom_col not found in file.")
    if(!start_col %in% names(cdict)) stop("start_col not found in file.")
    if(!end_col %in% names(cdict)) stop("end_col not found in file.")
    ### File MUST be bgzipped first
    messager("Ensuring file is bgzipped.",v=verbose)
    bgz_file <- Rsamtools::bgzip(file = path,
                                 dest = sprintf("%s.bgz",
                                                sub("\\.gz$|\\.bgz$", "",
                                                    path)),
                                 overwrite = overwrite)
    ### Tabix-index file
    messager("Tabix-indexing file.",v=verbose)  
    tbi_file <- Rsamtools::indexTabix(file = bgz_file, 
                                      seq = cdict[[chrom_col]],
                                      start = cdict[[start_col]],
                                      end = cdict[[end_col]],
                                      comment = names(cdict)[1]) 
    if(isTRUE(remove_tmp) && 
       file.exists(path)){
        messager("Removing temporary .tsv file.",v=verbose)  
        out <- file.remove(path)
    }
    return(bgz_file)
}
