# Check that FRQ column refers to minor/effect allele frequency not major

Check that FRQ column refers to minor/effect allele frequency not major

## Usage

``` r
check_frq_maf(sumstats_dt, frq_is_maf)
```

## Arguments

- frq_is_maf:

  Conventionally the FRQ column is intended to show the minor/effect
  allele frequency (MAF) but sometimes the major allele frequency can be
  inferred as the FRQ column. This logical variable indicates that the
  FRQ column should be renamed to MAJOR_ALLELE_FRQ if the frequency
  values appear to relate to the major allele i.e. \>0.5. By default
  this mapping won't occur i.e. is TRUE.

## Value

sumstats_dt, the modified summary statistics data table object
