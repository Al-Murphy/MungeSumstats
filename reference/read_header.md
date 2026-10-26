# Read in file header

Read in file header

## Usage

``` r
read_header(path, n = 2L, skip_vcf_metadata = FALSE, nThread = 1)
```

## Arguments

- path:

  Filepath for the summary statistics file to be formatted. A dataframe
  or datatable of the summary statistics file can also be passed
  directly to MungeSumstats using the path parameter.

- n:

  integer. The (maximal) number of lines to read. Negative values
  indicate that one should read up to the end of input on the
  connection.

- skip_vcf_metadata:

  logical, should VCF metadata be ignored

- nThread:

  Number of threads to use for parallel processes.

## Value

First `n` lines of the VCF header

## Examples

``` r
path <- system.file("extdata", "eduAttainOkbay.txt", 
                    package = "MungeSumstats") 
header <- read_header(path = path)                    
#> Reading header.
```
