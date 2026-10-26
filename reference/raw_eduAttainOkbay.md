# GWAS Educational Attainment Okbay 2016 - Subset

GWAS Summary Statistics on Educational Attainment by Okbay et al 2016:
PMID: 27898078 PMCID: PMC5509058 DOI: 10.1038/ng1216-1587b. A subset of
93 SNPs

## Format

txt document with 94 items

## Source

The summary statistics file was downloaded from
https://www.nature.com/articles/ng.3552 and formatted to a .rda with the
following:
` #Get example dataset, use Educational-Attainment_Okbay_2016 link<-"Educational-Attainment_Okbay_2016/EduYears_Discovery_5000.txt" eduAttainOkbay<-readLines(link,n=100) #There is an issue where values end with .0, this 0 is removed in func #There are also SNPs not on ref genome or arebi/tri allelic #So need to remove these in this dataset as its used for testing tmp <- tempfile() writeLines(eduAttainOkbay,con=tmp) eduAttainOkbay <- data.table::fread(tmp) #DT read removes the .0's #remove those not on ref genome and withbi/tri allelic rmv <- c("rs192818565","rs79925071","rs1606974","rs1871109", "rs73074378","rs7955289") eduAttainOkbay <- eduAttainOkbay[!MarkerName data.table::fwrite(eduAttainOkbay,file=tmp,sep="\t") eduAttainOkbay <- readLines(tmp) writeLines(eduAttainOkbay,"inst/extdata/eduAttainOkbay.txt") `

## Details

GWAS Summary Statistics on Educational Attainment by Okbay et al 2016
has been subsetted here to act as an example summary statistic file
which has some issues in the formatting. MungeSumstats can correct these
issues.

## eduAttainOkbay.txt

NA
