# Download VCF file and its index file from Open GWAS

Ideally, we would use [gwasvcf](https://github.com/MRCIEU/gwasvcf)
instead but it hasn't been made available on CRAN or Bioconductor yet,
so we can't include it as a dep.

## Usage

``` r
download_vcf(
  vcf_urls,
  vcf_dir = tempdir(),
  vcf_download = TRUE,
  download_method = "download.file",
  force_new = FALSE,
  quiet = FALSE,
  timeout = 10 * 60,
  nThread = 1
)
```

## Arguments

- vcf_urls:

  Remote URL to VCF file.

- vcf_dir:

  Where to download the original VCF from Open GWAS. *WARNING:* This is
  set to [`tempdir()`](https://rdrr.io/r/base/tempfile.html) by default.
  This means the raw (pre-formatted) VCFs be deleted upon ending the R
  session. Change this to keep the raw VCF file on disk (e.g.
  `vcf_dir="./raw_vcf"`).

- vcf_download:

  Download the original VCF from Open GWAS.

- download_method:

  `"axel"` (multi-threaded) or `"download.file"` (single-threaded) .

- force_new:

  Overwrite a previously downloaded VCF with the same path name.

- quiet:

  Run quietly.

- timeout:

  How many seconds before giving up on download. Passed to
  `download.file`. Default: `10*60` (10min).

- nThread:

  Number of threads to parallelize over.

## Value

List containing the paths to the downloaded VCF and its index file.

## Examples

``` r
#only run the examples if user has internet access:
if(try(is.character(getURL("www.google.com")))==TRUE){
vcf_urls <- ieugwasr::gwasinfo_files(id = "ieu-a-298")[[1]]
out_paths <- download_vcf(vcf_urls = vcf_urls)
}
#> Error in getURL("www.google.com") : could not find function "getURL"
```
