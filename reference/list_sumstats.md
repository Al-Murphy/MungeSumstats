# List munged summary statistics

Searches for and lists local GWAS summary statistics files munged by
[format_sumstats](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md)
or
[import_sumstats](https://al-murphy.github.io/MungeSumstats/reference/import_sumstats.md).

## Usage

``` r
list_sumstats(
  save_dir = getwd(),
  pattern = "*.tsv.gz$",
  ids_from_file = TRUE,
  verbose = TRUE
)
```

## Arguments

- save_dir:

  Top-level directory to recursively search for summary statistics files
  within.

- pattern:

  Regex pattern to search for files with.

- ids_from_file:

  Try to extract dataset IDs from file names. If `FALSE`, will infer IDs
  from the directory names instead.

- verbose:

  Print messages.

## Value

Named vector of summary stats paths.

## Examples

``` r
save_dir <- system.file("extdata",package = "MungeSumstats")
munged_files <- MungeSumstats::list_sumstats(save_dir = save_dir)
#> 1 file(s) found.
```
