# Convert summary statistics to desired object type

Convert summary statistics to desired object type

## Usage

``` r
convert_sumstats(
  sumstats_dt,
  return_format = c("data.table", "vranges", "granges")
)
```

## Arguments

- return_format:

  Object type to convert to; `"data.table"`, `"GenomicRanges"` or
  `"VRanges"`(default is `"data.table"`).

## Value

Summary statistics in the converted format
