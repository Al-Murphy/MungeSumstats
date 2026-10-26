# Ensure all SNPs have info score above threshold

Ensure all SNPs have info score above threshold

## Usage

``` r
check_info_score(
  sumstats_dt,
  INFO_filter,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files
)
```

## Arguments

- INFO_filter:

  numeric The minimum value permissible of the imputation information
  score (if present in sumstats file). Default 0.9.

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

  list of log file locations.

## Value

list containing sumstats_dt, the modified summary statistics data table
object and the log file list
