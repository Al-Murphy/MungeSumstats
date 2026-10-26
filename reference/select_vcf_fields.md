# Select VCF fields

Select non-empty columns from each VCF field type.

## Usage

``` r
select_vcf_fields(
  path,
  sampled_rows = 10000L,
  which = NULL,
  samples = NULL,
  nThread = 1,
  verbose = TRUE
)
```

## Arguments

- path:

  Path to local or remote VCF file.

- sampled_rows:

  First N rows to sample. Set `NULL` to use full `sumstats_file`. when
  determining whether cols are empty.

- which:

  Genomic ranges to be added if supplied. Default is NULL.

- samples:

  Which samples to use:

  - 1 : Only the first sample will be used (*DEFAULT*).

  - NULL : All samples will be used.

  - c("\<sample_id1\>","\<sample_id2\>",...) : Only user-selected
    samples will be used (case-insensitive).

- nThread:

  Number of threads to use for parallel processes.

- verbose:

  Print messages.

## Value

`ScanVcfParam` object.
