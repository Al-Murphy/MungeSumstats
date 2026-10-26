# Read VCF genome

Get the genome build of a remote or local VCF file.

## Usage

``` r
read_vcf_genome(
  header = NULL,
  validate = FALSE,
  default_genome = "HG19/GRCh37",
  verbose = TRUE
)
```

## Arguments

- header:

  Header extracted by
  [scanVcfHeader](https://rdrr.io/pkg/VariantAnnotation/man/scanVcf-methods.html).

- validate:

  Validate genome name using
  [mapGenomeBuilds](https://rdrr.io/pkg/GenomeInfoDb/man/mapGenomeBuilds.html).

- default_genome:

  When no genome can be extracted, default to this genome build.

- verbose:

  Print messages.

## Value

genome
