# GWAS Amyotrophic lateral sclerosis ieu open GWAS project - Subset

VCF (VCFv4.2) of the GWAS Amyotrophic lateral sclerosis ieu open GWAS
project Dataset: ebi-a-GCST005647. A subset of 99 SNPs

## Format

vcf document with 528 items relating to 99 SNPs

## Source

The summary statistics VCF (VCFv4.2) file was downloaded from
https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST005647/ and formatted to a
.rda with the following:
` #Get example VCF dataset, use GWAS Amyotrophic lateral sclerosis ALS_GWAS_VCF <- readLines("ebi-a-GCST005647.vcf.gz") #Subset to just the first 99 SNPs ALSvcf <- ALS_GWAS_VCF[1:528] writeLines(ALSvcf,"inst/extdata/ALSvcf.vcf") `

## Details

A VCF file (VCFv4.2) of the GWAS Amyotrophic lateral sclerosis ieu open
GWAS project has been subsetted here to act as an example summary
statistic file in VCF format which has some issues in the formatting.
MungeSumstats can correct these issues and produced a standardised
summary statistics format.

## ALSvcf.vcf

NA
