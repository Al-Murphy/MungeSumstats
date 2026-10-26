# axel downloader

R wrapper for axel, which enables multi-threaded download of a single
large file.

## Usage

``` r
axel(
  input_url,
  output_path,
  background = FALSE,
  nThread = 1,
  force_overwrite = FALSE,
  quiet = TRUE,
  alternate = TRUE,
  check_certificates = FALSE
)
```

## Arguments

- input_url:

  input_url.

- output_path:

  output_path.

- background:

  Run in background

- nThread:

  Number of threads to parallelize over.

- force_overwrite:

  Overwrite existing file.

- quiet:

  Run quietly.

- alternate:

  alternate,

- check_certificates:

  check_certificates

## Value

Path where the file has been downloaded

## See also

<https://github.com/axel-download-accelerator/axel/>

Other downloaders:
[`downloader()`](https://al-murphy.github.io/MungeSumstats/reference/downloader.md)
