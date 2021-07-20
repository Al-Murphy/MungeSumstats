#' Summary Statistics Column Headers
#'
#' @description List of uncorrected column headers often found in GWAS Summary 
#' Statistics column headers. Note the effect allele will always be the A2 
#' allele, this is the approach done for 
#' VCF(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7805039/). This is enforced
#' with the column header corrections here and also the check allele flipping
#' test.
#'
#' @source The code to prepare the .Rda file file from the marker file is:
#' \code{
#' # Most the data in the below table comes from the LDSC github wiki
#' data("sumstatsColHeaders")
#' # Make additions to sumstatsColHeaders using github version of MungeSumstats- 
#' # shown is an example of adding columns for Standard Error (SE)
#' #se_cols <- data.frame("Uncorrected"=c("SE","se","STANDARD.ERROR","STANDARD_ERROR","STANDARD-ERROR"), 
#' #                      "Corrected"=rep("SE",5))
#' #sumstatsColHeaders <- rbind(sumstatsColHeaders,se_cols)
#' #Once additions are made, save the new mapping dataset
#' usethis::use_data(sumstatsColHeaders,overwrite = TRUE, internal=TRUE)
#' save(sumstatsColHeaders,
#'       file="data/sumstatsColHeaders.rda")
#' # You will need to restart your r session for effects to take account      
#' }
#' @format dataframe with 99 rows and 2 columns
#' @usage data("sumstatsColHeaders")
"sumstatsColHeaders"


#' GWAS Educational Attainment Okbay 2016 - Subset
#'
#' @description GWAS Summary Statistics on Educational Attainment by Okbay et 
#' al 2016:
#' PMID: 27898078 PMCID: PMC5509058 DOI: 10.1038/ng1216-1587b.
#' A subset of 93 SNPs
#' 
#' @details GWAS Summary Statistics on Educational Attainment by Okbay et 
#' al 2016 has been subsetted here to act as an example summary statistic file
#' which has some issues in the formatting. MungeSumstats can correct these 
#' issues. 
#'
#' @name raw_eduAttainOkbay
#' @section eduAttainOkbay.txt
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
#' writeLines(eduAttainOkbay,"inst/extdata/eduAttainOkbay.txt")
#' }
#' @format txt document with 94 items
NULL


#' GWAS Amyotrophic lateral sclerosis ieu open GWAS project - Subset
#'
#' @description VCF (VCFv4.2) of the GWAS Amyotrophic lateral sclerosis ieu 
#' open GWAS project Dataset: ebi-a-GCST005647.
#' A subset of 99 SNPs
#'
#' @name raw_ALSvcf
#' @section ALSvcf.vcf
#' @details A VCF file (VCFv4.2) of the GWAS Amyotrophic lateral sclerosis ieu 
#' open GWAS project has been subsetted here to act as an example summary 
#' statistic file in VCF format which has some issues in the formatting. 
#' MungeSumstats can correct these issues and produced a standardised summary 
#' statistics format. 
#' @source The summary statistics VCF (VCFv4.2) file was downloaded from
#' https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST005647/
#' and formatted to a .rda with the following:
#' \code{
#' #Get example VCF dataset, use GWAS Amyotrophic lateral sclerosis
#' ALS_GWAS_VCF <- readLines("ebi-a-GCST005647.vcf.gz")
#' #Subset to just the first 99 SNPs
#' ALSvcf <- ALS_GWAS_VCF[1:528]
#' writeLines(ALSvcf,"inst/extdata/ALSvcf.vcf")
#' }
#' @format vcf document with 528 items relating to 99 SNPs
#' @details NULL
NULL

