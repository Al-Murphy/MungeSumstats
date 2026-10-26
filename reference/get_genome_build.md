# Infers the genome build of the summary statistics file (GRCh37 or GRCh38) from the data. Uses SNP (RSID) & CHR & BP to get genome build.

Infers the genome build of the summary statistics file (GRCh37 or
GRCh38) from the data. Uses SNP (RSID) & CHR & BP to get genome build.

## Usage

``` r
get_genome_build(
  sumstats,
  nThread = 1,
  sampled_snps = 10000,
  standardise_headers = TRUE,
  mapping_file = sumstatsColHeaders,
  dbSNP = 155,
  dbSNP_tarball = NULL,
  header_only = FALSE,
  allele_match_ref = FALSE,
  ref_genome = NULL,
  chr_filt = NULL
)
```

## Arguments

- sumstats:

  data table/data frame obj of the summary statistics file for the GWAS
  ,or file path to summary statistics file.

- nThread:

  Number of threads to use for parallel processes.

- sampled_snps:

  Downsample the number of SNPs used when inferring genome build to save
  time.

- standardise_headers:

  Run `standardise_sumstats_column_headers_crossplatform`.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in your file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See `data(sumstatsColHeaders)` for default mapping and
  necessary format.

- dbSNP:

  version of dbSNP to be used (144 or 155). Default is 155.

- dbSNP_tarball:

  Pass local versions of dbSNP in tarball format. Default of NULL uses
  the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball` was
  enabled to help with dbSNP versions \>=156, after the decision to no
  longer provide dbSNP releases as bioconductor packages. dbSNP 156
  tarball is available here: http://149.165.171.124/SNPlocs/.

- header_only:

  Instead of reading in the entire `sumstats` file, only read in the
  first N rows where N=`sampled_snps`. This should help speed up cases
  where you have to read in `sumstats` from disk each time.

- allele_match_ref:

  Instead of returning the genome_build this will return the proportion
  of matches to each genome build for each allele (A1,A2).

- ref_genome:

  name of the reference genome used for the GWAS ("GRCh37" or "GRCh38").
  Argument is case-insensitive. Default is NULL which infers the
  reference genome from the data.

- chr_filt:

  Internal for testing - filter reference genomes and sumstats to
  specific chromosomes for testing. Pass a list of chroms in format:
  c("1","2"). Default is NULL i.e. no filtering

## Value

ref_genome the genome build of the data
