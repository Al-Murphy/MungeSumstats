# Ensure that the first three columns are SNP, CHR, BP in that order and then A1, A2 if present

Ensure that the first three columns are SNP, CHR, BP in that order and
then A1, A2 if present

## Usage

``` r
check_col_order(sumstats_dt, path)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS

- path:

  Filepath for the summary statistics file to be formatted

## Value

list containing sumstats_dt, the modified summary statistics data table
object
