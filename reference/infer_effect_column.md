# Infer if effect relates to a1 or A2 if ambiguously named

Three checks are made to infer which allele the effect/frequency
information relates to if they are ambiguous (named A0, A1 and A2 or
equivalent):

1.  Check if ambiguous naming conventions are used (i.e. allele 0, 1 and
    2 or equivalent). If not exit, otherwise continue to next checks.
    This can be checked by using the mapping file and splitting A1/A2
    mappings by those that contain 0, 1 or 2 (ambiguous) or doesn't
    contain 0, 1 or 2 e.g. effect, tested (unambiguous so fine for MSS
    to handle as is).

2.  Look for effect column/frequency column where the A0/A1/A2
    explicitly mentioned, if found then we know the direction and should
    update A0/A1/A2 naming so A2 is the effect column. We can look for
    such columns by getting every combination of A0/A1/A2 naming and
    effect/frq naming.

3.  If not found in 2, a final check should be against the reference
    genome, whichever of A0, A1 and A2 has more of a match with the
    reference genome should be taken as **not** the effect allele. There
    is an assumption in this but is still better than guessing the
    ambiguous allele naming.

## Usage

``` r
infer_effect_column(
  sumstats_dt,
  dbSNP = 155,
  dbSNP_tarball = NULL,
  sampled_snps = 10000,
  mapping_file = sumstatsColHeaders,
  nThread = nThread,
  ref_genome = NULL,
  on_ref_genome = TRUE,
  infer_eff_direction = TRUE,
  eff_on_minor_alleles = FALSE,
  return_list = TRUE
)
```

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- dbSNP:

  version of dbSNP to be used for imputation (144 or 155). See
  `dbSNP_tarball` for different versions of dbSNP (including newer
  releases).

- dbSNP_tarball:

  Pass local versions of dbSNP in tarball format. Default of NULL uses
  the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball` was
  enabled to help with dbSNP versions \>=156, after the decision to no
  longer provide dbSNP releases as bioconductor packages. dbSNP 156
  tarball is available here: http://149.165.171.124/SNPlocs/.

- sampled_snps:

  Downsample the number of SNPs used when inferring genome build to save
  time.

- mapping_file:

  MungeSumstats has a pre-defined column-name mapping file which should
  cover the most common column headers and their interpretations.
  However, if a column header that is in youf file is missing of the
  mapping we give is incorrect you can supply your own mapping file.
  Must be a 2 column dataframe with column names "Uncorrected" and
  "Corrected". See data(sumstatsColHeaders) for default mapping and
  necessary format.

- nThread:

  Number of threads to use for parallel processes.

- ref_genome:

  name of the reference genome used for the GWAS ("GRCh37" or "GRCh38").
  Argument is case-insensitive. Default is NULL which infers the
  reference genome from the data.

- on_ref_genome:

  Binary Should a check take place that all SNPs are on the reference
  genome by SNP ID. Default is TRUE.

- infer_eff_direction:

  Binary Should a check take place to ensure the alleles match the
  effect direction? Default is TRUE.

- eff_on_minor_alleles:

  Binary Should MungeSumstats assume that the effects are majoritively
  measured on the minor alleles? Default is FALSE as this is an
  assumption that won't be appropriate in all cases. However, the
  benefit is that if we know the majority of SNPs have their effects
  based on the minor alleles, we can catch cases where the allele
  columns have been mislabelled.

- return_list:

  Return the `sumstats_dt` within a named list (default: `TRUE`).

## Value

list containing sumstats_dt, the modified summary statistics data table
object

## Details

Also, if eff_on_minor_alleles=TRUE, check 3 will be used in all cases.
However, This assumes that the effects are majoritively measured on the
minor alleles and should be used with caution as this is an assumption
that won't be appropriate in all cases. However, the benefit is that if
we know the majority of SNPs have their effects based on the minor
alleles, we can catch cases where the allele columns have been
mislabelled. IF eff_on_minor_alleles=TRUE, checks 1 and 2 will be
skipped.

## Examples

``` r
sumstats <- MungeSumstats::formatted_example()
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
#> Sorting coordinates with 'data.table'.
#for speed, don't run on_ref_genome part of check (on_ref_genome = FALSE)
sumstats_dt2<-infer_effect_column(sumstats_dt=sumstats,on_ref_genome = FALSE)
#> Infer Effect Column
#> First line of summary statistics file: 
#> SNP  CHR BP  A1  A2  FRQ BETA    SE  P   
#> Allele columns are ambiguous, attempting to infer direction
#> Can't infer allele columns from sumstats
```
