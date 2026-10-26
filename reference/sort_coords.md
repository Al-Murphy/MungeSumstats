# Sort sum stats

Sort summary statistics table by genomic coordinates.

## Usage

``` r
sort_coords(
  sumstats_dt,
  sort_coordinates = TRUE,
  sort_method = c("data.table", "GenomicRanges")
)
```

## Arguments

- sumstats_dt:

  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) obj
  of the summary statistics file for the GWAS.

- sort_method:

  Method to sort coordinates by:

  - "data.table" (default)Uses
    [setorderv](https://rdrr.io/pkg/data.table/man/setorder.html), which
    is must faster than "GenomicRanges" but less robust to variations in
    some sum stats files.

  - "GenomicRanges"Uses
    [sort.GenomicRanges](https://rdrr.io/pkg/GenomicRanges/man/GenomicRanges-comparison.html),
    which is more robust to variations in sum stats files but much
    slower than the "data.table" method.

- sort_coords:

  Whether to sort by coordinates.

## Value

Sorted sumstats_dt
