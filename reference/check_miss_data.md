# Remove SNPs with missing data

Remove SNPs with missing data

## Usage

``` r
check_miss_data(
  sumstats_dt,
  path,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files,
  drop_na_cols
)
```

## Arguments

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

- drop_na_cols:

  A character vector of column names to be checked for missing values.
  Rows with missing values in any of these columns (if present in the
  dataset) will be dropped. If `NULL`, all columns will be checked for
  missing values. Default columns are SNP, chromosome, position, allele
  1, allele2, effect columns (frequency, beta, Z-score, standard error,
  log odds, signed sumstats, odds ratio), p value and N columns.

## Value

list containing sumstats_dt, the modified summary statistics data table
object and a log file list.
