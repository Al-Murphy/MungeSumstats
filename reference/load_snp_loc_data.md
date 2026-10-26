# Loads the SNP locations and alleles for Homo sapiens from dbSNP builds

Loads the SNP locations and alleles for Homo sapiens from dbSNP builds

## Usage

``` r
load_snp_loc_data(ref_genome, dbSNP, dbSNP_tarball = NULL, msg = NULL)
```

## Arguments

- ref_genome:

  character, "GRCh37" or "GRCh38"

- dbSNP:

  integer, dbSNP build number (144, 155, or any installed SNPlocs
  package)

- dbSNP_tarball:

  Optional path to a .tar.gz containing: one or more .rds files (Bioc
  SNPlocs package layout).

- msg:

  optional character to message before loading

## Value

A data.table or OnDiskLongTable of SNP locations
