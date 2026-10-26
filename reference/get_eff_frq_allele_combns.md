# Get combinations of uncorrected allele and effect (and frq) columns

Get combinations of uncorrected allele and effect (and frq) columns

## Usage

``` r
get_eff_frq_allele_combns(
  mapping_file = sumstatsColHeaders,
  eff_frq_cols = c("BETA", "OR", "LOG_ODDS", "SIGNED_SUMSTAT", "Z", "FRQ")
)
```

## Arguments

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in youf file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

- eff_frq_cols:

  Corrected effect or frequency column names found in a sumstats.
  Default of BETA, OR, LOG_ODDS, SIGNED_SUMSTAT, Z and FRQ.

## Value

datatable containing uncorrected and corrected combinations
