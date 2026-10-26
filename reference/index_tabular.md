# Tabix-index file: table

Convert summary stats file to tabix format.

## Usage

``` r
index_tabular(
  path,
  chrom_col = "CHR",
  start_col = "BP",
  end_col = start_col,
  overwrite = TRUE,
  remove_tmp = TRUE,
  verbose = TRUE
)
```

## Source

Borrowed function from
[echotabix](https://github.com/RajLabMSSM/echotabix/blob/main/R/convert.R).

## Arguments

- path:

  Path to GWAS summary statistics file.

- chrom_col:

  Name of the chromosome column in `sumstats_dt` (e.g. "CHR").

- start_col:

  Name of the starting genomic position column in `sumstats_dt` (e.g.
  "POS","start").

- end_col:

  Name of the ending genomic position column in `sumstats_dt` (e.g.
  "POS","end"). Can be the same as `start_col` when `sumstats_dt` only
  contains SNPs that span 1 base pair (bp) each.

- overwrite:

  A logical(1) indicating whether `dest` should be over-written, if it
  already exists.

- remove_tmp:

  Remove the temporary uncompressed version of the file (*.tsv*).

- verbose:

  Print messages.

## Value

Path to tabix-indexed tabular file

## See also

Other tabix:
[`index_vcf()`](https://al-murphy.github.io/MungeSumstats/reference/index_vcf.md)

## Examples

``` r
sumstats_dt <- MungeSumstats::formatted_example() 
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
#> Sorting coordinates with 'data.table'.
path <- tempfile(fileext = ".tsv")
MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt, save_path = path)
#> Writing in tabular format ==> /tmp/Rtmpt871Bb/file3252740186a3.tsv
indexed_file <- MungeSumstats::index_tabular(path = path)
#> Converting full summary stats file to tabix format for fast querying...
#> Reading header.
#> Ensuring file is bgzipped.
#> Tabix-indexing file.
#> Removing temporary .tsv file.
```
