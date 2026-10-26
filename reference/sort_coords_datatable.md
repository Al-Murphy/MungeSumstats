# Sort sum stats: data.table

Sort summary statistics table by genomic coordinates using a fast
`data.table`-native strategy

## Usage

``` r
sort_coords_datatable(
  sumstats_dt,
  chr_col = "CHR",
  start_col = "BP",
  end_col = start_col
)
```

## Arguments

- sumstats_dt:

  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) obj
  of the summary statistics file for the GWAS.

- chr_col:

  Chromosome column name.

- start_col:

  Genomic end position column name.

## Value

Sorted sumstats_dt
