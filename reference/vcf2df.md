# VCF to DF

Function to convert a VariantAnnotation `CollapsedVCF`/`ExpandedVCF`
object to a `data.frame`.

## Usage

``` r
vcf2df(
  vcf,
  add_sample_names = TRUE,
  add_rowranges = TRUE,
  drop_empty_cols = TRUE,
  unique_cols = TRUE,
  unique_rows = TRUE,
  unlist_cols = TRUE,
  sampled_rows = NULL,
  verbose = TRUE
)
```

## Source

[Original code
source](https://gist.github.com/zhujack/849b75f5a8305edaeca1001dfb9c3fe9)

### vcfR

if(!require("pinfsc50")) install.packages("pinfsc50") vcf_file \<-
system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50") vcf \<-
read.vcfR( vcf_file, verbose = FALSE ) vcf_df_list \<-
vcfR::vcfR2tidy(vcf, single_frame=TRUE) vcf_df \<-
data.table::data.table(vcf_df_list\$dat)

## Arguments

- vcf:

  Variant Call Format (VCF) file imported into R as a VariantAnnotation
  [CollapsedVCF](https://rdrr.io/pkg/VariantAnnotation/man/VCF-class.html)/
  [ExpandedVCF](https://rdrr.io/pkg/VariantAnnotation/man/VCF-class.html)
  object.

- add_sample_names:

  Append sample names to column names (e.g. "EZ" –\> "EZ_ubm-a-2929").

- add_rowranges:

  Include `rowRanges` from VCF as well.

- drop_empty_cols:

  Drop columns that are filled entirely with: `NA`, `"."`, or `""`.

- unique_cols:

  Only keep uniquely named columns.

- unique_rows:

  Only keep unique rows.

- unlist_cols:

  If any columns are lists instead of vectors, unlist them. Required to
  be `TRUE` when `unique_rows=TRUE`.

- sampled_rows:

  First N rows to sample. Set `NULL` to use full `sumstats_file`. when
  determining whether cols are empty.

- verbose:

  Print messages.

## Value

data.frame version of VCF

## Examples

``` r
  
#### VariantAnnotation ####
# path <- "https://github.com/brentp/vcfanno/raw/master/example/exac.vcf.gz"
path <- system.file("extdata", "ALSvcf.vcf",
                    package = "MungeSumstats")
                    
vcf <- VariantAnnotation::readVcf(file = path)
vcf_df <- MungeSumstats:::vcf2df(vcf = vcf)
#> Converting VCF to data.table.
#> Expanding VCF first, so number of rows may increase.
#> Checking for empty columns.
#> Removing 2 empty columns.
#> Unlisting 4 columns.
#> Dropped 314 duplicate rows.
#> Time difference of 0.1 secs
#> VCF data.table contains: 101 rows x 12 columns.
```
