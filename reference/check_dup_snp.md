# Ensure all rows have unique SNP IDs, drop those that don't

Ensure all rows have unique SNP IDs, drop those that don't

## Usage

``` r
check_dup_snp(
  sumstats_dt,
  indels,
  path,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files,
  bi_allelic_filter,
  check_dups
)
```

## Arguments

- indels:

  Binary does your Sumstats file contain Indels? These don't exist in
  our reference file so they will be excluded from checks if this value
  is TRUE. Default is TRUE.

- path:

  Filepath for the summary statistics file to be formatted. A dataframe
  or datatable of the summary statistics file can also be passed
  directly to MungeSumstats using the path parameter.

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

- bi_allelic_filter:

  Binary Should non-bi-allelic SNPs be removed. Default is TRUE.

- check_dups:

  whether to check for duplicates - if formatting QTL datasets this
  should be set to FALSE otherwise keep as TRUE. Default is TRUE.

## Value

list containing sumstats_dt, the modified summary statistics data table
object and log files list
