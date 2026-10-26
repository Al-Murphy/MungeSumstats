# Map column names to positions.

Useful in situations where you need to specify columns by index instead
of name (e.g. awk queries).

## Usage

``` r
column_dictionary(file_path)
```

## Source

Borrowed function from
[echotabix](https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R).

` eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats" ) tmp <- tempfile(fileext = ".tsv") file.copy(eduAttainOkbayPth, tmp) cdict <- MungeSumstats:::column_dictionary(file_path = tmp) `

## Arguments

- file_path:

  Path to full summary stats file (or any really file you want to make a
  column dictionary for).

## Value

Named list of column positions.
