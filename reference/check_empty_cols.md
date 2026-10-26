# Check for empty columns

Empty columns contain only ".", NA, or 0

## Usage

``` r
check_empty_cols(sumstats_dt, sampled_rows = NULL, verbose = TRUE)
```

## Arguments

- sampled_rows:

  First N rows to sample. Set `NULL` to use full `sumstats_file`. when
  determining whether cols are empty.

- verbose:

  Print messages.

## Value

empty_cols
