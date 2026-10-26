# Ensure that the N column is all integers

Ensure that the N column is all integers

## Usage

``` r
check_n_int(sumstats_dt, path, convert_n_int, imputation_ind)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS

- path:

  Filepath for the summary statistics file to be formatted

- convert_n_int:

  Binary, if N (the number of samples) is not an integer, should this be
  rounded? Default is TRUE.

- imputation_ind:

  Binary Should a column be added for each imputation step to show what
  SNPs have imputed values for differing fields. This includes a field
  denoting SNP allele flipping (flipped). **Note** these columns will be
  in the formatted summary statistics returned. Default is FALSE.

## Value

list containing sumstats_dt, the modified summary statistics data table
object.
