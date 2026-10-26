# Remove non-biallelic SNPs

Remove non-biallelic SNPs

## Usage

``` r
check_bi_allelic(
  sumstats_dt,
  path,
  ref_genome,
  bi_allelic_filter,
  rsids,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files,
  dbSNP,
  dbSNP_tarball
)
```

## Arguments

- path:

  Filepath for the summary statistics file to be formatted. A dataframe
  or datatable of the summary statistics file can also be passed
  directly to MungeSumstats using the path parameter.

- ref_genome:

  name of the reference genome used for the GWAS ("GRCh37" or "GRCh38").
  Argument is case-insensitive. Default is NULL which infers the
  reference genome from the data.

- bi_allelic_filter:

  Binary Should non-bi-allelic SNPs be removed. Default is TRUE.

- log_folder_ind:

  Binary Should log files be stored containing all filtered out SNPs
  (separate file per filter). The data is outputted in the same format
  specified for the resulting sumstats file. The only exception to this
  rule is if output is vcf, then log file saved as .tsv.gz. Default is
  FALSE.

- tabix_index:

  Index the formatted summary statistics with
  [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.

- nThread:

  Number of threads to use for parallel processes.

- log_files:

  list of log file locations

- dbSNP:

  version of dbSNP to be used for imputation (144 or 155). See
  `dbSNP_tarball` for different versions of dbSNP (including newer
  releases).

- dbSNP_tarball:

  Pass local versions of dbSNP in tarball format. Default of NULL uses
  the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball` was
  enabled to help with dbSNP versions \>=156, after the decision to no
  longer provide dbSNP releases as bioconductor packages. dbSNP 156
  tarball is available here: http://149.165.171.124/SNPlocs/.

## Value

A list containing two data tables:

- `sumstats_dt`: the modified summary statistics data table object

- `rsids`: snpsById, filtered to SNPs of interest if loaded already. Or
  else `NULL`.

- `log_files`: log file list
