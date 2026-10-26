# To `GRanges`

Convert a
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html) to
[GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html).

## Usage

``` r
to_granges(
  sumstats_dt,
  seqnames.field = "CHR",
  start.field = "BP",
  end.field = "BP",
  style = c("NCBI", "UCSC")
)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- seqnames.field:

  A character vector of recognized names for the column in `df` that
  contains the chromosome name (a.k.a. sequence name) associated with
  each genomic range. Only the first name in `seqnames.field` that is
  found in `colnames(df)` is used. If no one is found, then an error is
  raised.

- start.field:

  A character vector of recognized names for the column in `df` that
  contains the start positions of the genomic ranges. Only the first
  name in `start.field` that is found in `colnames(df)` is used. If no
  one is found, then an error is raised.

- end.field:

  A character vector of recognized names for the column in `df` that
  contains the end positions of the genomic ranges. Only the first name
  in `start.field` that is found in `colnames(df)` is used. If no one is
  found, then an error is raised.

- style:

  `GRanges` style to convert to, "NCBI" or "UCSC".

## Value

`GRanges` object
