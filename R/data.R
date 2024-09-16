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
#' # Shown is an example of adding new A1 and A2 naming
#' a1_name <- c("NON","RISK","ALLELE")
#' a2_name <- c("RISK","ALLELE")
#' all_delims <- c("_",".",""," ","-")
#' all_uncorr_a1 <- vector(mode="list",length = length(all_delims))
#' all_corr_a1 <- vector(mode="list",length = length(all_delims))
#' all_uncorr_a2 <- vector(mode="list",length = length(all_delims))
#' all_corr_a2 <- vector(mode="list",length = length(all_delims))
#' for(i in seq_along(all_delims)){
#' delim <- all_delims[i]
#' a1 <- unlist(paste(a1_name,collapse=delim))
#' a2 <- unlist(paste(a2_name,collapse=delim))
#' all_uncorr_a1[[i]] <- a1
#' all_uncorr_a2[[i]] <- a2
#' all_corr_a1[[i]] <- "A1"
#'   all_corr_a2[[i]] <- "A2"
#' }
#' se_cols <- data.frame("Uncorrected"=c(unlist(all_uncorr_a1),unlist(all_uncorr_a2)),
#'                      "Corrected"=c(unlist(all_corr_a1),unlist(all_corr_a2)))
#' # Or another example .....
#' # shown is an example of adding columns for Standard Error (SE)
#' se_cols <- data.frame("Uncorrected"=c("SE","se","STANDARD.ERROR",
#'                                       "STANDARD_ERROR","STANDARD-ERROR"),
#'                      "Corrected"=rep("SE",5))
#' sumstatsColHeaders <- rbind(sumstatsColHeaders,se_cols)
#' #Once additions are made, order & save the new mapping dataset
#' #now sort ordering -important for logic that 
#' # uncorrected=corrected comes first
#' sumstatsColHeaders$ordering <-
#'     sumstatsColHeaders$Uncorrected==sumstatsColHeaders$Corrected
#' sumstatsColHeaders <-
#'     sumstatsColHeaders[order(sumstatsColHeaders$Corrected,
#'                              sumstatsColHeaders$ordering,decreasing = TRUE),]
#' rownames(sumstatsColHeaders)<-1:nrow(sumstatsColHeaders)
#' sumstatsColHeaders$ordering <- NULL
#' #manually move FREQUENCY to above MAR - github issue 95
#' frequency <- sumstatsColHeaders[sumstatsColHeaders$Uncorrected=="FREQUENCY",]
#' maf <- sumstatsColHeaders[sumstatsColHeaders$Uncorrected=="MAF",]
#' if(as.integer(rownames(frequency))>as.integer(rownames(maf))){
#'   sumstatsColHeaders[as.integer(rownames(frequency)),] <- maf
#'   sumstatsColHeaders[as.integer(rownames(maf)),] <- frequency
#' }   
#' usethis::use_data(sumstatsColHeaders,overwrite = TRUE, internal=TRUE)
#' save(sumstatsColHeaders,
#'       file="data/sumstatsColHeaders.rda")
#' # You will need to restart your r session for effects to take account
#' }
#' @format dataframe with 2 columns
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

#' UCSC Chain file hg38 to hg19
#'
#' @description UCSC Chain file hg38 to hg19, .chain.gz file, downloaded from 
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/ on 09/10/21
#'
#' @name hg38ToHg19
#' @section hg38ToHg19.over.chain.gz
#' @details UCSC Chain file hg38 to hg19, .chain.gz file, downloaded on 09/10/21
#' To be used as a back up if the download from UCSC fails.
#' @source The chain file was downloaded from
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/ 
#' \code{
#' utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',tempdir())
#' }
#' @format gunzipped chain file
#' @details NULL
NULL

#' UCSC Chain file hg19 to hg38
#'
#' @description UCSC Chain file hg19 to hg38, .chain.gz file, downloaded from 
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/ on 09/10/21
#'
#' @name hg19ToHg38
#' @section hg19ToHg38.over.chain.gz
#' @details UCSC Chain file hg19 to hg38, .chain.gz file, downloaded on 09/10/21
#' To be used as a back up if the download from UCSC fails.
#' @source The chain file was downloaded from
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/ 
#' \code{
#' utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',tempdir())
#' }
#' @format gunzipped chain file
#' @details NULL
NULL

#' Local ieu-a-298 file from IEU Open GWAS
#'
#' @description Local ieu-a-298 file from IEU Open GWAS, downloaded on 09/10/21.
#'
#' @name ieu-a-298
#' @section ieu-a-298.tsv.gz
#' @details Local ieu-a-298 file from IEU Open GWAS, downlaoded on 09/10/21. 
#' This is done in case the download in the package vignette fails.
#' @source The file was downloaded with:
#' \code{
#' MungeSumstats::import_sumstats(ids = "ieu-a-298",ref_genome = "GRCH37")
#' }
#' @format gunzipped tsv file
#' @details NULL
NULL
