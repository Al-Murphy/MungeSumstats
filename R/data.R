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

#' GWAS Educational Attainment Okbay 2016 - Subset
#'
#' GWAS Summary Statistics on Educational Attainment by Okbay et al 2016:
#' PMID: 27898078 PMCID: PMC5509058 DOI: 10.1038/ng1216-1587b.
#' A subset of 99 SNPs
#'
#' @source The summary statistics file was downloaded from
#' https://www.nature.com/articles/ng.3552
#' and formatted to a .rda with the following:
#' \code{
#' #Get example dataset, use Educational-Attainment_Okbay_2016
#' link<-"Educational-Attainment_Okbay_2016/EduYears_Discovery_5000.txt"
#' eduAttainOkbay<-readLines(link,n=100)
#' #There is an issue where values end with .0, this 0 is removed in func
#' #So need to remove it in this dataset as its used for testing
#' tmp <- tempfile()
#' writeLines(eduAttainOkbay,con=tmp)
#' eduAttainOkbay <- data.table::fread(tmp) #DT read removes the .0's
#' data.table::fwrite(eduAttainOkbay,file=tmp,sep="\t")
#' eduAttainOkbay <- readLines(tmp)
#' usethis::use_data(eduAttainOkbay,overwrite = TRUE)
#' }
#' @format character vector with 100 items
"eduAttainOkbay"
