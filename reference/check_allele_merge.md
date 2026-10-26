# Ensure that A1:A2 or A1/A2 or A1\>A2 or A2\>A1 aren't merged into 1 column

Ensure that A1:A2 or A1/A2 or A1\>A2 or A2\>A1 aren't merged into 1
column

## Usage

``` r
check_allele_merge(sumstats_dt, path)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS

- path:

  Filepath for the summary statistics file to be formatted

## Value

list containing sumstats_dt, the modified summary statistics data table
object.
