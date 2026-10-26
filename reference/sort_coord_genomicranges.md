# Sort sum stats: GenomicRanges

Sort summary statistics table by genomic coordinates using a slower (but
in some cases more robust) `GenomicRanges` strategy

## Usage

``` r
sort_coord_genomicranges(sumstats_dt)
```

## Arguments

- sumstats_dt:

  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) obj
  of the summary statistics file for the GWAS.

## Value

Sorted sumstats_dt
