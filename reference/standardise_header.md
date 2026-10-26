# Standardise the column headers in the Summary Statistics files

Use a reference data table of common column header names (stored in
`sumstatsColHeaders` or user inputted mapping file) to convert them to a
standard set, i.e. chromosome -\> CHR. This function does not check that
all the required column headers are present. The amended header is
written directly back into the file

## Usage

``` r
standardise_header(
  sumstats_dt,
  mapping_file = sumstatsColHeaders,
  uppercase_unmapped = TRUE,
  convert_A0 = TRUE,
  return_list = TRUE
)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in youf file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

- uppercase_unmapped:

  For columns that could not be identified in the `mapping_file`, return
  them in the same format they were input as (without forcing them to
  uppercase).

- convert_A0:

  Whether to convert A\* (representing A0) to A1/A2. This should be done
  unless checking if A0 was present in the input as if you do it you
  can't infer this. Default is TRUE

- return_list:

  Return the `sumstats_dt` within a named list (default: `TRUE`).

## Value

list containing sumstats_dt, the modified summary statistics data table
object

## Examples

``` r
sumstats_dt <- data.table::fread(system.file("extdata", "eduAttainOkbay.txt",
                                             package = "MungeSumstats"))
sumstats_dt2 <- standardise_header(sumstats_dt=sumstats_dt)
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
```
