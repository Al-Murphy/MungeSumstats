#' Determine summary statistics file type and read them into memory
#'
#' @return \code{data.table} of formatted summary statistics
#'
#' @param nrows integer. The (maximal) number of lines to read.
#' If \code{Inf}, will read in all rows.
#' @param standardise_headers Standardise headers first.
#' @inheritParams format_sumstats
#' @inheritParams standardise_sumstats_column_headers_crossplatform
#' @export
#' @importFrom data.table fread
#' @examples
#' path <- system.file("extdata", "eduAttainOkbay.txt",
#'     package = "MungeSumstats"
#' )
#' eduAttainOkbay <- read_sumstats(path = path)
read_sumstats <- function(path,
                          nThread = 1,
                          nrows = Inf,
                          standardise_headers = FALSE,
                          mapping_file = sumstatsColHeaders) {
    if (is.data.frame(path)) {
        message("Summary statistics passed as R object.")
        sumstats_file <- data.table::as.data.table(path)
        if (!is.infinite(nrows)) {
            sumstats_file <- sumstats_file[seq(1, nrows), ]
        }
    } else {
        header <- read_header(path = path)
        is_vcf <- check_vcf(header = header)
        if (is_vcf) {
            message("Importing VCF file: ", path)
            sumstats_file <- read_vcf(path = path, nThread = nThread)
        } else {
            is_tabular <- check_tabular(header = header)
            if (is_tabular) {
                message("Importing tabular file: ", path)
                sumstats_file <- data.table::fread(path,
                    nThread = nThread,
                    nrows = nrows
                )
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
    if (standardise_headers) {
        sumstats_file <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_file,
                mapping_file = mapping_file
            )[["sumstats_dt"]]
    }
    return(sumstats_file)
}
