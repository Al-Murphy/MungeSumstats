#' Map column names to positions.
#'
#' Useful in situations where you need to specify columns by
#' index instead of name (e.g. awk queries).
#'
#' @source Borrowed function from
#' \href{https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R}{
#' echotabix}.
#' 
#' @source 
#' \code{
#' eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
#'     package = "MungeSumstats"
#' )
#' tmp <- tempfile(fileext = ".tsv")
#' file.copy(eduAttainOkbayPth, tmp)
#' cdict <- MungeSumstats:::column_dictionary(file_path = tmp)
#' }
#' 
#' @param file_path Path to full summary stats file
#' (or any really file you want to make a column dictionary for).
#' 
#' @return Named list of column positions.
#' 
#' @keywords internal
#' @importFrom stats setNames
column_dictionary <- function(file_path) {
    # Get the index of each column name
    header <- read_header(path = file_path, 
                          # n must be 2 or else 
                          # fread won't be able to parse text
                          n = 2, 
                          skip_vcf_metadata = TRUE) 
    cNames <- colnames(data.table::fread(text = header))
    colDict <- stats::setNames(seq(1, length(cNames)), cNames)
    return(colDict)
}
