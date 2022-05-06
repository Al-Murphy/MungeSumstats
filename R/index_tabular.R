#' Tabix-index file: table
#' 
#' Convert summary stats file to tabix format.
#'
#' @source Borrowed function from 
#' \href{https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R}{
#' echotabix}.
#' 
#' @param path Path to GWAS summary statistics file.
#' @param verbose Print messages.
#' @inheritParams liftover
#' 
#' @return Path to tabix-indexed tabular file
#'
#' @family tabix
#' @examples
#' sumstats_dt <- MungeSumstats::formatted_example()
#' sumstats_dt <- MungeSumstats:::sort_coords(sumstats_dt = sumstats_dt)
#' path <- tempfile(fileext = ".tsv")
#' MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt, save_path = path)
#' 
#' indexed_file <- MungeSumstats::index_tabular(path = path)
#' @export
#' @importFrom R.utils gunzip isGzipped
index_tabular <- function(path,
                          chrom_col = "CHR",
                          start_col = "BP",
                          end_col = start_col,
                          verbose = TRUE) {
    msg <- paste("Converting full summary stats file to",
                 "tabix format for fast querying...")
    message(msg)
    #### Read header and make dictionary ####
    cdict <- column_dictionary(file_path = path)
    #### Check column exist ####
    if(!chrom_col %in% names(cdict)) stop("chrom_col not found in file.")
    if(!start_col %in% names(cdict)) stop("start_col not found in file.")
    if(!end_col %in% names(cdict)) stop("end_col not found in file.")
    #### Make sure input file isn't empty ####
    if (file.size(path) == 0) {
        msg2 <- paste("Removing empty file:", path)
        message(msg2)
        out <- file.remove(path)
    }
    #### Ensure it's not already compressed ####
    if(R.utils::isGzipped(path)){
        path <- R.utils::gunzip(path, overwrite = TRUE)
    }
    ### File MUST be bgzipped first
    message("Ensuring file is bgzipped.")
    bgz_file <- Rsamtools::bgzip(file = path,
                                 dest = sprintf("%s.bgz",
                                                sub("\\.gz$|\\.bgz$", "",
                                                    path)),
                                 overwrite = TRUE)
    ### Tabix-index file
    message("Tabix-indexing file.") 
    # Rsamtools::indexTabix is not user-friendly, 
    # and is prone to errors without additional wrappers functions
    #(ie seqminer::tabix.createIndex)
    # tbx <- Rsamtools::indexTabix(file = bgz_file,
    #                               seq = chrom_col,
    #                               start = start_col,
    #                               comment = "SNP", 
    #                               end = end_col)
    seqminer::tabix.createIndex(
        bgzipFile = bgz_file,
        sequenceColumn = cdict[[chrom_col]],
        startColumn = cdict[[start_col]],
        endColumn = cdict[[end_col]],
        ## Just use the first columns name
        metaChar = names(cdict)[1]
    )
    return(bgz_file)
}