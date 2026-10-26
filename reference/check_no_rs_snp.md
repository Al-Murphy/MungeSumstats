# Ensure that SNP appears to be valid RSIDs (starts with rs)

Ensure that SNP appears to be valid RSIDs (starts with rs)

## Usage

``` r
check_no_rs_snp(
  sumstats_dt,
  path,
  ref_genome,
  snp_ids_are_rs_ids,
  indels,
  imputation_ind,
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

- snp_ids_are_rs_ids:

  Binary Should the supplied SNP ID's be assumed to be RSIDs. If not,
  imputation using the SNP ID for other columns like base-pair position
  or chromosome will not be possible. If set to FALSE, the SNP RS ID
  will be imputed from the reference genome if possible. Default is
  TRUE.

- indels:

  Binary does your Sumstats file contain Indels? These don't exist in
  our reference file so they will be excluded from checks if this value
  is TRUE. Default is TRUE.

- imputation_ind:

  Binary Should a column be added for each imputation step to show what
  SNPs have imputed values for differing fields. This includes a field
  denoting SNP allele flipping (flipped). On the flipped value, this
  denoted whether the alelles where switched based on MungeSumstats
  initial choice of A1, A2 from the input column headers and thus may
  not align with what the creator intended.**Note** these columns will
  be in the formatted summary statistics returned. Default is FALSE.

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

list containing sumstats_dt, the modified summary statistics data table
object and the log file list.
