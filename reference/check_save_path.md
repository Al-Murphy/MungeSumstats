# Check if save path and log folder is appropriate

Check if save path and log folder is appropriate

## Usage

``` r
check_save_path(
  save_path,
  log_folder,
  log_folder_ind,
  tabix_index,
  write_vcf = FALSE,
  verbose = TRUE
)
```

## Arguments

- save_path:

  File path to save formatted data. Defaults to
  `tempfile(fileext=".tsv.gz")`.

- log_folder:

  Filepath to the directory for the log files and the log of
  MungeSumstats messages to be stored. Default is a temporary directory.
  Note the name of the log files (log messages and log outputs) are now
  the same as the name of the file specified in the save path parameter
  with the extension '\_log_msg.txt' and '\_log_output.txt'
  respectively.

- log_folder_ind:

  Binary Should log files be stored containing all filtered out SNPs
  (separate file per filter). The data is outputted in the same format
  specified for the resulting sumstats file. The only exception to this
  rule is if output is vcf, then log file saved as .tsv.gz. Default is
  FALSE.

- tabix_index:

  Index the formatted summary statistics with
  [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.

- write_vcf:

  Whether to write as VCF (TRUE) or tabular file (FALSE).

- verbose:

  Print messages.

## Value

Corrected `save_path`, the file type, the separator, corrected
`log_folder`,the log file extension.
