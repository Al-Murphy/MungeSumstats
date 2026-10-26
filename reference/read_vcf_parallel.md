# Read VCF: parallel

Read a VCF file across 1 or more threads in parallel. If `tilewidth` is
not specified, the size of each chunk will be determined by total genome
size divided by `ntile`. By default, `ntile` is equal to the number of
threads, `nThread`. For further discussion on how this function was
optimised, see
[here](https://github.com/Bioconductor/VariantAnnotation/issues/59) and
[here](https://github.com/Bioconductor/VariantAnnotation/issues/57).

## Usage

``` r
read_vcf_parallel(
  path,
  samples = 1,
  which = NULL,
  use_params = TRUE,
  as_datatable = TRUE,
  sampled_rows = 10000L,
  include_xy = FALSE,
  download = TRUE,
  vcf_dir = tempdir(),
  download_method = "download.file",
  force_new = FALSE,
  tilewidth = NULL,
  mt_thresh = 100000L,
  nThread = 1,
  ntile = nThread,
  verbose = TRUE
)
```

## Source

` path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz" #### Single-threaded #### vcf <- MungeSumstats:::read_vcf_parallel(path = path) #### Parallel #### vcf2 <- MungeSumstats:::read_vcf_parallel(path = path, nThread=11) `

## Arguments

- path:

  Path to local or remote VCF file.

- samples:

  Which samples to use:

  - 1 : Only the first sample will be used (*DEFAULT*).

  - NULL : All samples will be used.

  - c("\<sample_id1\>","\<sample_id2\>",...) : Only user-selected
    samples will be used (case-insensitive).

- which:

  Genomic ranges to be added if supplied. Default is NULL.

- use_params:

  When `TRUE` (default), increases the speed of reading in the VCF by
  omitting columns that are empty based on the head of the VCF (NAs
  only). NOTE that that this requires the VCF to be sorted,
  bgzip-compressed, tabix-indexed, which
  [read_vcf](https://al-murphy.github.io/MungeSumstats/reference/read_vcf.md)
  will attempt to do.

- as_datatable:

  Return the data as a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  (default: `TRUE`) or a
  [VCF](https://rdrr.io/pkg/VariantAnnotation/man/VCF-class.html)
  (`FALSE`).

- sampled_rows:

  First N rows to sample. Set `NULL` to use full `sumstats_file`. when
  determining whether cols are empty.

- download:

  Download the VCF (and its index file) to a temp folder before reading
  it into R. This is important to keep `TRUE` when `nThread>1` to avoid
  making too many queries to remote file.

- vcf_dir:

  Where to download the original VCF from Open GWAS. *WARNING:* This is
  set to [`tempdir()`](https://rdrr.io/r/base/tempfile.html) by default.
  This means the raw (pre-formatted) VCFs be deleted upon ending the R
  session. Change this to keep the raw VCF file on disk (e.g.
  `vcf_dir="./raw_vcf"`).

- download_method:

  `"axel"` (multi-threaded) or `"download.file"` (single-threaded) .

- force_new:

  If a formatted file of the same names as `save_path` exists,
  formatting will be skipped and this file will be imported instead
  (default). Set `force_new=TRUE` to override this.

- tilewidth:

  The desired tile width. The effective tile width might be slightly
  different but is guaranteed to never be more than the desired width.

- mt_thresh:

  When the number of rows (variants) in the VCF is `< mt_thresh`, only
  use single-threading for reading in the VCF. This is because the
  overhead of parallelisation outweighs the speed benefits when VCFs are
  small.

- nThread:

  Number of threads to use for parallel processes.

- ntile:

  The number of tiles to generate.

- verbose:

  Print messages.

## Value

VCF file.
