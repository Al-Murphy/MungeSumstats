# Read -log10 p-value column

Parse p-value column in VCF file.of other general -loq10 p-values

## Usage

``` r
read_log_pval(
  sumstats_dt,
  mapping_file = sumstatsColHeaders,
  return_list = TRUE
)
```

## Arguments

- sumstats_dt:

  Summary stats data.table.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in youf file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

- return_list:

  Binary, whether to return the dt in a list or not - list is standard
  for the
  [`format_sumstats()`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md)
  function.

## Value

Null output.
