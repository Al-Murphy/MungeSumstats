# Example logs file

Example logs file produced by
[format_sumstats](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md).

## Usage

``` r
logs_example(read = FALSE)
```

## Source

` eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats") sumstats_dt <- data.table::fread(eduAttainOkbayPth) #### Introduce values that need to be fixed #### sumstats_dt$Pval[10:15] <- 5 sumstats_dt$Pval[20:22] <- -5 sumstats_dt$Pval[23:25] <- "5e-324" ss_path <- tempfile() data.table::fwrite(sumstats_dt, ss_path) log_folder <- tempdir() reformatted <- MungeSumstats::format_sumstats( path = ss_path, ref_genome = "GRCh37", log_folder = log_folder, log_mungesumstats_msgs = TRUE, log_folder_ind = TRUE, ) file.copy(reformatted$log_files$MungeSumstats_log_msg, "inst/extdata",overwrite = TRUE) `

## Arguments

- read:

  Whether to read the logs file into memory.

## Value

Path to logs file.
