# Check for N column if not present and user wants, impute N based on user's sample size. **NOTE** this will be the same value for each SNP which is not necessarily correct and may cause issues down the line. N can also be inputted with "ldsc", "sum", "giant" or "metal" by passing one or multiple of these.

Check for N column if not present and user wants, impute N based on
user's sample size. **NOTE** this will be the same value for each SNP
which is not necessarily correct and may cause issues down the line. N
can also be inputted with "ldsc", "sum", "giant" or "metal" by passing
one or multiple of these.

## Usage

``` r
compute_nsize(
  sumstats_dt,
  imputation_ind = FALSE,
  compute_n = c("ldsc", "giant", "metal", "sum"),
  standardise_headers = FALSE,
  force_new = FALSE,
  return_list = TRUE
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

- compute_n:

  How to compute per-SNP sample size (new column "N").

  - `0`: N will not be computed.

  - `>0`: If any number \>0 is provided, that value will be set as N for
    every row. **Note**: Computing N this way is incorrect and should be
    avoided if at all possible.

  - `"sum"`: N will be computed as: cases (N_CAS) + controls (N_CON), so
    long as both columns are present.

  - `"ldsc"`: N will be computed as effective sample size: Neff
    =(N_CAS+N_CON)\*(N_CAS/(N_CAS+N_CON)) /
    mean((N_CAS/(N_CAS+N_CON))(N_CAS+N_CON)==max(N_CAS+N_CON)).

  - `"giant"`: N will be computed as effective sample size: Neff = 2 /
    (1/N_CAS + 1/N_CON).

  - `"metal"`: N will be computed as effective sample size: Neff = 4 /
    (1/N_CAS + 1/N_CON).

- standardise_headers:

  Standardise headers first.

- force_new:

  If "Neff" (or "N") already exists in `sumstats_dt`, replace it with
  the recomputed version.

- return_list:

  Return the `sumstats_dt` within a named list (default: `TRUE`).

## Value

`list("sumstats_dt"=sumstats_dt)`

## Examples

``` r
sumstats_dt <- MungeSumstats::formatted_example()
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
#> Sorting coordinates with 'data.table'.
sumstats_dt2 <- MungeSumstats::compute_nsize(sumstats_dt=sumstats_dt,
                                             compute_n=10000)
#> Assigning N=10000 for all SNPs.
```
