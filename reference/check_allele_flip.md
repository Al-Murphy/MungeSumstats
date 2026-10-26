# Ensure A1 & A2 are correctly named, if GWAS SNP constructed as Alternative/Reference or Risk/Nonrisk alleles these SNPs will need to be converted to Reference/Alternative or Nonrisk/Risk. Here non-risk is defined as what's on the reference genome (this may not always be the case).

Ensure A1 & A2 are correctly named, if GWAS SNP constructed as
Alternative/Reference or Risk/Nonrisk alleles these SNPs will need to be
converted to Reference/Alternative or Nonrisk/Risk. Here non-risk is
defined as what's on the reference genome (this may not always be the
case).

## Usage

``` r
check_allele_flip(
  sumstats_dt,
  path,
  ref_genome,
  rsids,
  allele_flip_check,
  allele_flip_drop,
  allele_flip_z,
  allele_flip_frq,
  bi_allelic_filter,
  flip_frq_as_biallelic,
  imputation_ind,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files,
  standardise_headers = FALSE,
  mapping_file,
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

- allele_flip_check:

  Binary Should the allele columns be checked against reference genome
  to infer if flipping is necessary. Default is TRUE.

- allele_flip_drop:

  Binary Should the SNPs for which neither their A1 or A2 base pair
  values match a reference genome be dropped. Default is TRUE.

- allele_flip_z:

  Binary should the Z-score be flipped along with effect and FRQ columns
  like Beta? It is assumed to be calculated off the effect size not the
  P-value and so will be flipped i.e. default TRUE.

- allele_flip_frq:

  Binary should the frequency (FRQ) column be flipped along with effect
  and z-score columns like Beta? Default TRUE.

- bi_allelic_filter:

  Binary Should non-bi-allelic SNPs be removed. Default is TRUE.

- flip_frq_as_biallelic:

  Binary Should non-bi-allelic SNPs frequency values be flipped as 1-p
  despite there being other alternative alleles? Default is FALSE but if
  set to TRUE, this allows non-bi-allelic SNPs to be kept despite
  needing flipping.

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

- standardise_headers:

  Run `standardise_sumstats_column_headers_crossplatform` first.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in youf file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

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

- `sumstats_dt`: the modified summary statistics `data.table` object.

- `rsids`: snpsById, filtered to SNPs of interest if loaded already. Or
  else NULL.

- `log_files`: log file list
