# Ensure that there is at least one signed column in summary statistics file Impute beta if user requests

Ensure that there is at least one signed column in summary statistics
file Impute beta if user requests

## Usage

``` r
check_signed_col(
  sumstats_dt,
  impute_beta,
  log_folder_ind,
  rsids,
  imputation_ind,
  check_save_out,
  tabix_index,
  log_files,
  nThread
)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS

- impute_beta:

  Binary, whether BETA should be imputed using other effect data if it
  isn't present in the sumstats. Note that this imputation is an
  approximation (for Z & SE approach) so could have an effect on
  downstream analysis. Use with caution. The different methods
  MungeSumstats will try and impute beta (in this order or priority)
  are:

  1.  log(OR) 2. Z x SE Default value is FALSE.

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

- log_files:

  list of log file locations

- nThread:

  Number of threads to use for parallel processes.

## Value

null
