# Check numeric columns

Checks for any columns that should be numeric, and ensures that they are
indeed numeric.

## Usage

``` r
check_numeric(sumstats_dt, cols = c("P", "SE", "FRQ", "MAF", "BETA"))
```

## Arguments

- sumstats_dt:

  Summary stats with column names already standardised by
  [format_sumstats](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md).

- cols:

  Names of columns that should be numeric. If any of these columns are
  not actually present in `sumstats_dt`, they will be skipped.

## Value

sumstats_dt
