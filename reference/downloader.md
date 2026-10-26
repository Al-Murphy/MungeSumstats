# downloader wrapper

R wrapper for [axel](https://github.com/axel-download-accelerator/axel)
(multi-threaded) and
[download.file](https://rdrr.io/r/utils/download.file.html)
(single-threaded) download functions.

## Usage

``` r
downloader(
  input_url,
  output_path,
  download_method = "axel",
  background = FALSE,
  force_overwrite = FALSE,
  quiet = TRUE,
  show_progress = TRUE,
  continue = TRUE,
  nThread = 1,
  alternate = TRUE,
  check_certificates = TRUE,
  timeout = 10 * 60
)
```

## Source

[Suggestion to avoid 'proc\$get_built_file() : Build process
failed'](https://stackoverflow.com/a/66602496/13214824)

## Arguments

- input_url:

  input_url.

- output_path:

  output_path.

- download_method:

  `"axel"` (multi-threaded) or `"download.file"` (single-threaded) .

- background:

  Run in background

- force_overwrite:

  Overwrite existing file.

- quiet:

  Run quietly.

- show_progress:

  show_progress.

- continue:

  continue.

- nThread:

  Number of threads to parallelize over.

- alternate:

  alternate,

- check_certificates:

  check_certificates

- timeout:

  How many seconds before giving up on download. Passed to
  `download.file`. Default: `10*60` (10min).

## Value

Local path to downloaded file.

## See also

Other downloaders:
[`axel()`](https://al-murphy.github.io/MungeSumstats/reference/axel.md)
