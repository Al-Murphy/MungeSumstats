# Register cores

Register a multi-threaded instances using BiocParallel.

## Usage

``` r
register_cores(workers = 1, progressbar = TRUE)
```

## Arguments

- workers:

  `integer(1)` Number of workers. Defaults to the maximum of 1 or the
  number of cores determined by `detectCores` minus 2 unless environment
  variables `R_PARALLELLY_AVAILABLECORES_FALLBACK` or
  `BIOCPARALLEL_WORKER_NUMBER` are set otherwise. For a `SOCK` cluster,
  `workers` can be a
  [`character()`](https://rdrr.io/r/base/character.html) vector of host
  names.

- progressbar:

  `logical(1)` Enable progress bar (based on plyr:::progress_text).

## Value

Null output.
