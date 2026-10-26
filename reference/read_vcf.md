# Read in VCF file

Read in a VCF file as a
[VCF](https://rdrr.io/pkg/VariantAnnotation/man/VCF-class.html) or a
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html). Can
optionally save the VCF/data.table as well.

## Usage

``` r
read_vcf(
  path,
  as_datatable = TRUE,
  save_path = NULL,
  tabix_index = FALSE,
  samples = 1,
  which = NULL,
  use_params = TRUE,
  sampled_rows = 10000L,
  download = TRUE,
  vcf_dir = tempdir(),
  download_method = "download.file",
  force_new = FALSE,
  mt_thresh = 100000L,
  nThread = 1,
  verbose = TRUE
)
```

## Source

` `[`#### Benchmarking #### library(VCFWrenchR) library(VariantAnnotation) path <- "https://gwas.mrcieu.ac.uk/files/ubm-a-2929/ubm-a-2929.vcf.gz" vcf <- VariantAnnotation::readVcf(file = path) N <- 1e5 vcf_sub <- vcf[1:N,] res <- microbenchmark::microbenchmark( "vcf2df"={dat1 <- MungeSumstats:::vcf2df(vcf = vcf_sub)}, "VCFWrenchR"= {dat2 <- as.data.frame(x = vcf_sub)}, "VRanges"={dat3 <- data.table::as.data.table( methods::as(vcf_sub, "VRanges"))}, times=1 )`](https://rdrr.io/r/base/library.html)` `

[Discussion on VariantAnnotation
GitHub](https://github.com/Bioconductor/VariantAnnotation/issues/57)

[Discussion on VariantAnnotation
GitHub](https://github.com/Bioconductor/VariantAnnotation/issues/59)

## Arguments

- path:

  Path to local or remote VCF file.

- as_datatable:

  Return the data as a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  (default: `TRUE`) or a
  [VCF](https://rdrr.io/pkg/VariantAnnotation/man/VCF-class.html)
  (`FALSE`).

- save_path:

  File path to save formatted data. Defaults to
  `tempfile(fileext=".tsv.gz")`.

- tabix_index:

  Index the formatted summary statistics with
  [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.

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
  bgzip-compressed, tabix-indexed, which read_vcf will attempt to do.

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

- mt_thresh:

  When the number of rows (variants) in the VCF is `< mt_thresh`, only
  use single-threading for reading in the VCF. This is because the
  overhead of parallelisation outweighs the speed benefits when VCFs are
  small.

- nThread:

  Number of threads to use for parallel processes.

- verbose:

  Print messages.

## Value

The VCF file in data.table format.

## Examples

``` r
#### Local file ####
path <- system.file("extdata","ALSvcf.vcf", package="MungeSumstats")
sumstats_dt <- read_vcf(path = path)
#> Loading required namespace: GenomicFiles
#> Using local VCF.
#> bgzip-compressing VCF file.
#> Finding empty VCF columns based on first 10,000 rows.
#> Dropping 1 duplicate column(s).
#> 1 sample detected: EBI-a-GCST005647
#> Constructing ScanVcfParam object.
#> VCF contains: 39,630,630 variant(s) x 1 sample(s)
#> Reading VCF file: single-threaded
#> Converting VCF to data.table.
#> Expanding VCF first, so number of rows may increase.
#> Dropping 1 duplicate column(s).
#> Checking for empty columns.
#> Unlisting 3 columns.
#> Dropped 314 duplicate rows.
#> Time difference of 0 secs
#> VCF data.table contains: 101 rows x 11 columns.
#> Time difference of 0.2 secs
#> Renaming ID as SNP.
#> sumstats has -log10 P-values; these will be converted to unadjusted p-values in the 'P' column.
#> No INFO (SI) column detected.

#### Remote file ####
## Small GWAS (0.2Mb)
# path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz"
# sumstats_dt2 <- read_vcf(path = path)

## Large GWAS (250Mb)
# path <- "https://gwas.mrcieu.ac.uk/files/ubm-a-2929/ubm-a-2929.vcf.gz"
# sumstats_dt3 <- read_vcf(path = path, nThread=11)

### Very large GWAS (500Mb)
# path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-1124/ieu-a-1124.vcf.gz"
# sumstats_dt4 <- read_vcf(path = path, nThread=11)
```
