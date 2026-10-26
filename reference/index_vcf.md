# Tabix-index file: VCF

Convert summary stats file to tabix format

## Usage

``` r
index_vcf(path, verbose = TRUE)
```

## Source

Borrowed function from
[echotabix](https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R).

## Arguments

- path:

  Path to VCF.

- verbose:

  Print messages.

## Value

Path to tabix-indexed tabular file

## See also

Other tabix:
[`index_tabular()`](https://al-murphy.github.io/MungeSumstats/reference/index_tabular.md)

## Examples

``` r
eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
                                 package = "MungeSumstats")
sumstats_dt <- data.table::fread(eduAttainOkbayPth, nThread = 1)
sumstats_dt <- 
MungeSumstats:::standardise_sumstats_column_headers_crossplatform(
    sumstats_dt = sumstats_dt)$sumstats_dt
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
sumstats_dt <- MungeSumstats:::sort_coords(sumstats_dt = sumstats_dt)
#> Sorting coordinates with 'data.table'.
path <- tempfile(fileext = ".tsv")
MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt, save_path = path)
#> Writing in tabular format ==> /tmp/Rtmpt871Bb/file3252391a3db.tsv
    
indexed_file <- MungeSumstats::index_tabular(path = path)
#> Converting full summary stats file to tabix format for fast querying...
#> Reading header.
#> Ensuring file is bgzipped.
#> Tabix-indexing file.
#> Removing temporary .tsv file.
```
