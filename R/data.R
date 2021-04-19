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
#' A subset of 93 SNPs
#'
#' @source The summary statistics file was downloaded from
#' https://www.nature.com/articles/ng.3552
#' and formatted to a .rda with the following:
#' \code{
#' #Get example dataset, use Educational-Attainment_Okbay_2016
#' link<-"Educational-Attainment_Okbay_2016/EduYears_Discovery_5000.txt"
#' eduAttainOkbay<-readLines(link,n=100)
#' #There is an issue where values end with .0, this 0 is removed in func
#' #There are also SNPs not on ref genome or arebi/tri allelic
#' #So need to remove these in this dataset as its used for testing
#' tmp <- tempfile()
#' writeLines(eduAttainOkbay,con=tmp)
#' eduAttainOkbay <- data.table::fread(tmp) #DT read removes the .0's
#' #remove those not on ref genome and withbi/tri allelic
#' rmv <- c("rs192818565","rs79925071","rs1606974","rs1871109",
#'          "rs73074378","rs7955289")
#' eduAttainOkbay <- eduAttainOkbay[!MarkerName %in% rmv,]
#' data.table::fwrite(eduAttainOkbay,file=tmp,sep="\t")
#' eduAttainOkbay <- readLines(tmp)
#' usethis::use_data(eduAttainOkbay,overwrite = TRUE)
#' }
#' @format character vector with 94 items
"eduAttainOkbay"

#' GWAS Amyotrophic lateral sclerosis ieu open GWAS project - Subset
#'
#' VCF (VCFv4.2) of the GWAS Amyotrophic lateral sclerosis ieu open GWAS project
#' Dataset: ebi-a-GCST005647.
#' A subset of 99 SNPs
#'
#' @source The summary statistics VCF (VCFv4.2) file was downloaded from
#' https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST005647/
#' and formatted to a .rda with the following:
#' \code{
#' #Get example VCF dataset, use GWAS Amyotrophic lateral sclerosis
#' AML_GWAS_VCF <- readLines("ebi-a-GCST005647.vcf.gz")
#' #Subset to just the first 99 SNPs
#' ieuAmlVcf <- AML_GWAS_VCF[1:528]
#' usethis::use_data(ieuAmlVcf,overwrite = TRUE)
#' }
#' @format character vector with 528 items relating to 99 SNPs
"ieuAmlVcf"
