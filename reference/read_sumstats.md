# Determine summary statistics file type and read them into memory

Determine summary statistics file type and read them into memory

## Usage

``` r
read_sumstats(
  path,
  nrows = Inf,
  standardise_headers = FALSE,
  samples = 1,
  sampled_rows = 10000L,
  nThread = 1,
  mapping_file = sumstatsColHeaders
)
```

## Arguments

- path:

  Filepath for the summary statistics file to be formatted. A dataframe
  or datatable of the summary statistics file can also be passed
  directly to MungeSumstats using the path parameter.

- nrows:

  integer. The (maximal) number of lines to read. If `Inf`, will read in
  all rows.

- standardise_headers:

  Standardise headers first.

- samples:

  Which samples to use:

  - 1 : Only the first sample will be used (*DEFAULT*).

  - NULL : All samples will be used.

  - c("\<sample_id1\>","\<sample_id2\>",...) : Only user-selected
    samples will be used (case-insensitive).

- sampled_rows:

  First N rows to sample. Set `NULL` to use full `sumstats_file`. when
  determining whether cols are empty.

- nThread:

  Number of threads to use for parallel processes.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in youf file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

## Value

`data.table` of formatted summary statistics

## Examples

``` r
path <- system.file("extdata", "eduAttainOkbay.txt",
    package = "MungeSumstats"
)
eduAttainOkbay <- read_sumstats(path = path)
#> Importing tabular file: /home/runner/work/_temp/Library/MungeSumstats/extdata/eduAttainOkbay.txt
#> Checking for empty columns.
```
