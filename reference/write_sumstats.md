# Write sum stats file to disk

Write sum stats file to disk

## Usage

``` r
write_sumstats(
  sumstats_dt,
  save_path,
  ref_genome = NULL,
  sep = "\t",
  write_vcf = FALSE,
  save_format = NULL,
  tabix_index = FALSE,
  nThread = 1,
  return_path = FALSE,
  save_path_check = FALSE
)
```

## Source

[VariantAnnotation::writeVcf has some unexpected/silent file renaming
behavior](https://github.com/Bioconductor/VariantAnnotation/issues/35)

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- save_path:

  File path to save formatted data. Defaults to
  `tempfile(fileext=".tsv.gz")`.

- ref_genome:

  name of the reference genome used for the GWAS ("GRCh37" or "GRCh38").
  Argument is case-insensitive. Default is NULL which infers the
  reference genome from the data.

- sep:

  The separator between columns. Defaults to the character in the set
  `[,\t |;:]` that separates the sample of rows into the most number of
  lines with the same number of fields. Use `NULL` or `""` to specify no
  separator; i.e. each line a single character column like
  [`base::readLines`](https://rdrr.io/r/base/readLines.html) does.

- write_vcf:

  Whether to write as VCF (TRUE) or tabular file (FALSE).

- save_format:

  Output format of sumstats. Options are NULL - standardised output
  format from MungeSumstats, LDSC - output format compatible with LDSC
  and openGWAS - output compatible with openGWAS VCFs. Default is NULL.
  **NOTE** - If LDSC format is used, the naming convention of A1 as the
  reference (genome build) allele and A2 as the effect allele will be
  reversed to match LDSC (A1 will now be the effect allele). See more
  info on this
  [here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68). Note
  that any effect columns (e.g. Z) will be inrelation to A1 now instead
  of A2.

- tabix_index:

  Index the formatted summary statistics with
  [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.

- nThread:

  The number of threads to use. Experiment to see what works best for
  your data on your hardware.

- return_path:

  Return `save_path`. This will have been modified in some cases (e.g.
  after compressing and tabix-indexing a previously un-compressed file).

- save_path_check:

  Ensure path name is valid (given the other arguments) before writing
  (default: FALSE).

## Value

If `return_path=TRUE`, returns `save_path`. Else returns `NULL`.

## Examples

``` r
path <- system.file("extdata", "eduAttainOkbay.txt",
    package = "MungeSumstats"
)
eduAttainOkbay <- read_sumstats(path = path)
#> Importing tabular file: /home/runner/work/_temp/Library/MungeSumstats/extdata/eduAttainOkbay.txt
#> Checking for empty columns.
write_sumstats(
    sumstats_dt = eduAttainOkbay,
    save_path = tempfile(fileext = ".tsv.gz")
)
#> Writing in tabular format ==> /tmp/Rtmpt871Bb/file32527a2b547c.tsv.gz
```
