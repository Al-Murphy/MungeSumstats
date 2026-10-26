# Ensure all rows are unique based on SNP,CHR,BP,A1,A2, drop those that aren't

Ensure all rows are unique based on SNP,CHR,BP,A1,A2, drop those that
aren't

## Usage

``` r
check_dup_row(
  sumstats_dt,
  check_dups,
  path,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files
)
```

## Arguments

- check_dups:

  whether to check for duplicates - if formatting QTL datasets this
  should be set to FALSE otherwise keep as TRUE. Default is TRUE.

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

## Value

list containing sumstats_dt, the modified summary statistics data table
object and log files list
