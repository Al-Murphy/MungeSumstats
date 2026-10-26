# Ensure that the standard error (se) is positive for all SNPs Also impute se if missing

Ensure that the standard error (se) is positive for all SNPs Also impute
se if missing

## Usage

``` r
check_pos_se(
  sumstats_dt,
  path,
  pos_se,
  log_folder_ind,
  imputation_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files,
  impute_se
)
```

## Arguments

- path:

  Filepath for the summary statistics file to be formatted. A dataframe
  or datatable of the summary statistics file can also be passed
  directly to MungeSumstats using the path parameter.

- pos_se:

  Binary Should the standard Error (SE) column be checked to ensure it
  is greater than 0? Those that are, are removed (if present in sumstats
  file). Default TRUE.

- log_folder_ind:

  Binary Should log files be stored containing all filtered out SNPs
  (separate file per filter). The data is outputted in the same format
  specified for the resulting sumstats file. The only exception to this
  rule is if output is vcf, then log file saved as .tsv.gz. Default is
  FALSE.

- imputation_ind:

  Binary Should a column be added for each imputation step to show what
  SNPs have imputed values for differing fields. This includes a field
  denoting SNP allele flipping (flipped). On the flipped value, this
  denoted whether the alelles where switched based on MungeSumstats
  initial choice of A1, A2 from the input column headers and thus may
  not align with what the creator intended.**Note** these columns will
  be in the formatted summary statistics returned. Default is FALSE.

- tabix_index:

  Index the formatted summary statistics with
  [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.

- nThread:

  Number of threads to use for parallel processes.

- log_files:

  list of log file locations

- impute_se:

  Binary, whether the standard error should be imputed using other
  effect data if it isn't present in the sumstats. Note that this
  imputation is an approximation so could have an effect on downstream
  analysis. Use with caution. The different methods MungeSumstats will
  try and impute se (in this order or priority) are:

  1.  BETA / Z 2. abs(BETA/ qnorm(P/2)) Default is FALSE.

## Value

list containing sumstats_dt, the modified summary statistics data table
object and the log file list
