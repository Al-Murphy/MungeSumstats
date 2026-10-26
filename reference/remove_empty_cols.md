# Remove empty columns

Remote columns that are empty or contain all the same values in a
data.table.

## Usage

``` r
remove_empty_cols(sumstats_dt, sampled_rows = NULL, verbose = TRUE)
```

## Arguments

- sampled_rows:

  First N rows to sample. Set `NULL` to use full `sumstats_file`. when
  determining whether cols are empty.

- verbose:

  Print messages.

## Value

Null output.
