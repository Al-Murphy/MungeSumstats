# Formatted example

Returns an example of summary stats that have had their column names
already standardised with
[standardise_header](https://al-murphy.github.io/MungeSumstats/reference/standardise_header.md).

## Usage

``` r
formatted_example(
  path = system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats"),
  formatted = TRUE,
  sorted = TRUE
)
```

## Arguments

- path:

  Path to raw example file. Default to built-in dataset.

- formatted:

  Whether the column names should be formatted (default:`TRUE`).

- sorted:

  Whether the rows should be sorted by genomic coordinates
  (default:`TRUE`).

## Value

`sumstats_dt`

## Examples

``` r
sumstats_dt <- MungeSumstats::formatted_example()
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
#> Sorting coordinates with 'data.table'.
```
