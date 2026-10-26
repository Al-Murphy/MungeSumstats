# Ensure that the standard error (se) is positive for all SNPs

Ensure that the standard error (se) is positive for all SNPs

## Usage

``` r
check_effect_columns_nonzero(
  sumstats_dt,
  path,
  effect_columns_nonzero,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files
)
```

## Arguments

- path:

  Filepath for the summary statistics file to be formatted. A dataframe
  or datatable of the summary statistics file can also be passed
  directly to MungeSumstats using the path parameter.

- effect_columns_nonzero:

  Binary should the effect columns in the data BETA,OR (odds
  ratio),LOG_ODDS,SIGNED_SUMSTAT be checked to ensure no SNP=0. Those
  that do are removed(if present in sumstats file). Default FALSE.

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
object and the log file list
