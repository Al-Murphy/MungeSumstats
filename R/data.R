#' Summary Statistics Column Headers
#'
#' List of uncorrected column headers often found in GWAS Summary Statistics 
#' column headers
#'
#' @source The code to prepare the .Rda file file from the marker file is:
#' \code{
#' # Most the data in the below table comes from the LDSC github wiki
#' sumstatsColHeaders <- read.csv("inst/extdata/Magma_Column_headers.csv",
#'      stringsAsFactors = FALSE)
#' usethis::use_data(sumstatsColHeaders,overwrite = TRUE)
#' }
#' @format dataframe with 82 rows nd 2 columns
"sumstatsColHeaders"