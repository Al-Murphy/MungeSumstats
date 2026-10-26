# Check for Z-score column

The following ensures that a Z-score column is present. The Z-score
formula we used here is a R implementation of the formula used in
[LDSC's
munge_sumstats.py](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L363):

## Usage

``` r
check_zscore(
  sumstats_dt,
  imputation_ind,
  compute_z = "BETA",
  force_new_z = FALSE,
  standardise_headers = FALSE,
  mapping_file
)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- imputation_ind:

  Binary Should a column be added for each imputation step to show what
  SNPs have imputed values for differing fields. This includes a field
  denoting SNP allele flipping (flipped). **Note** these columns will be
  in the formatted summary statistics returned. Default is FALSE.

- compute_z:

  Whether to compute Z-score column. Default is FALSE. This can be
  computed from Beta and SE with (Beta/SE) or P
  (Z:=sign(BETA)\*sqrt(stats::qchisq(P,1,lower=FALSE))). **Note** that
  imputing the Z-score from P for every SNP will not be perfectly
  correct and may result in a loss of power. This should only be done as
  a last resort. Use 'BETA' to impute by BETA/SE and 'P' to impute by
  SNP p-value.

- force_new_z:

  When a "Z" column already exists, it will be used by default. To
  override and compute a new Z-score column from P set
  `force_new_z=TRUE`.

- standardise_headers:

  Run `standardise_sumstats_column_headers_crossplatform` first.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in your file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

## Value

`list("sumstats_dt"=sumstats_dt)`

## Details

`np.sqrt(chi2.isf(P, 1))`

The R implementation is adapted from the `GenomicSEM::munge` function,
after optimizing for speed using `data.table`:

`sumstats_dt[,Z:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))]`

*NOTE*: `compute_z` is set to `TRUE` by default to ensure
standardisation of the "Z" column (which can be computed differently in
different datasets).
