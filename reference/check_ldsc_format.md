# Ensures that parameters are compatible with LDSC format

Format summary statistics for direct input to Linkage Disequilibrium
SCore (LDSC) regression without the need to use their
`munge_sumstats.py` script first.

## Usage

``` r
check_ldsc_format(
  sumstats_dt,
  save_format,
  convert_n_int,
  allele_flip_check,
  compute_z,
  compute_n
)
```

## Source

[LDSC GitHub](https://github.com/bulik/ldsc)

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- save_format:

  Output format of sumstats. Options are NULL - standardised output
  format from MungeSumstats, LDSC - output format compatible with LDSC
  and openGWAS - output compatible with openGWAS VCFs. Default is NULL.
  **NOTE** - If LDSC format is used, the naming convention of A1 as the
  reference (genome build) allele and A2 as the effect allele will be
  reversed to match LDSC (A1 will now be the effect allele). See more
  info on this
  [here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68). Note
  that any effect columns (e.g. Z) will be inrelation to A1 now instead
  of A2.

- convert_n_int:

  Binary, if N (the number of samples) is not an integer, should this be
  rounded? Default is TRUE.

- allele_flip_check:

  Binary Should the allele columns be checked against reference genome
  to infer if flipping is necessary. Default is TRUE.

- compute_z:

  Whether to compute Z-score column. Default is FALSE. This can be
  computed from Beta and SE with (Beta/SE) or P
  (Z:=sign(BETA)\*sqrt(stats::qchisq(P,1,lower=FALSE))). **Note** that
  imputing the Z-score from P for every SNP will not be perfectly
  correct and may result in a loss of power. This should only be done as
  a last resort. Use 'BETA' to impute by BETA/SE and 'P' to impute by
  SNP p-value.

- compute_n:

  Whether to impute N. Default of 0 won't impute, any other integer will
  be imputed as the N (sample size) for every SNP in the dataset.
  **Note** that imputing the sample size for every SNP is not correct
  and should only be done as a last resort. N can also be inputted with
  "ldsc", "sum", "giant" or "metal" by passing one of these for this
  field or a vector of multiple. Sum and an integer value creates an N
  column in the output whereas giant, metal or ldsc create an Neff or
  effective sample size. If multiples are passed, the formula used to
  derive it will be indicated.

## Value

Formatted summary statistics

## Details

[LDSC
documentation](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format).
