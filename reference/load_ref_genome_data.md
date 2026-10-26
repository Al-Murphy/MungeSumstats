# Load the reference genome data for SNPs of interest

Load the reference genome data for SNPs of interest

## Usage

``` r
load_ref_genome_data(
  snps,
  ref_genome,
  dbSNP = c(144, 155),
  dbSNP_tarball = NULL,
  msg = NULL,
  chr_filt = NULL
)
```

## Source

` sumstats_dt <- formatted_example() rsids <- MungeSumstats:::load_ref_genome_data(snps = sumstats_dt$SNP, ref_genome = "GRCH37", dbSNP=144) `

## Arguments

- snps:

  Character vector SNPs by rs_id from sumstats file of interest.

- ref_genome:

  Name of the reference genome used for the GWAS (GRCh37 or GRCh38)

- dbSNP:

  version of dbSNP to be used (144 or 155)

- dbSNP_tarball:

  Pass local versions of dbSNP in tarball format. Default of NULL uses
  the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball` was
  enabled to help with dbSNP versions \>=156, after the decision to no
  longer provide dbSNP releases as bioconductor packages. dbSNP 156
  tarball is available here: http://149.165.171.124/SNPlocs/.

- msg:

  Optional name of the column missing from the dataset in question.
  Default is NULL

- chr_filt:

  Internal for testing - filter reference genomes and sumstats to
  specific chromosomes for testing. Pass a list of chroms in format:
  c("1","2"). Default is NULL i.e. no filtering.

## Value

data table of snpsById, filtered to SNPs of interest.
