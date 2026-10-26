# Genome build liftover

Transfer genomic coordinates from one genome build to another.

## Usage

``` r
liftover(
  sumstats_dt,
  convert_ref_genome,
  ref_genome,
  chain_source = "ensembl",
  imputation_ind = TRUE,
  chrom_col = "CHR",
  start_col = "BP",
  end_col = start_col,
  as_granges = FALSE,
  style = "NCBI",
  local_chain = NULL,
  rmv_chr = c(),
  verbose = TRUE
)
```

## Source

[liftOver](https://doi.org/doi:10.18129/B9.bioc.liftOver)

[UCSC chain
files](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/)

[Ensembl chain
files](https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/)

## Arguments

- sumstats_dt:

  data table obj of the summary statistics file for the GWAS.

- convert_ref_genome:

  name of the reference genome to convert to ("GRCh37" or "GRCh38").
  This will only occur if the current genome build does not match.
  Default is not to convert the genome build (NULL).

- ref_genome:

  name of the reference genome used for the GWAS ("GRCh37" or "GRCh38").
  Argument is case-insensitive. Default is NULL which infers the
  reference genome from the data.

- chain_source:

  chain file source used ("ucsc" as default, or "ensembl")

- imputation_ind:

  Binary Should a column be added for each imputation step to show what
  SNPs have imputed values for differing fields. This includes a field
  denoting SNP allele flipping (flipped). On the flipped value, this
  denoted whether the alelles where switched based on MungeSumstats
  initial choice of A1, A2 from the input column headers and thus may
  not align with what the creator intended.**Note** these columns will
  be in the formatted summary statistics returned. Default is FALSE.

- chrom_col:

  Name of the chromosome column in `sumstats_dt` (e.g. "CHR").

- start_col:

  Name of the starting genomic position column in `sumstats_dt` (e.g.
  "POS","start").

- end_col:

  Name of the ending genomic position column in `sumstats_dt` (e.g.
  "POS","end"). Can be the same as `start_col` when `sumstats_dt` only
  contains SNPs that span 1 base pair (bp) each.

- as_granges:

  Return results as
  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  instead of a
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  (default: `FALSE`).

- style:

  Style to return
  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object in (e.g. "NCBI" = 4; "UCSC" = "chr4";) (default: `"NCBI"`).

- local_chain:

  Path to local chain file to use instead of downlaoding. Default of
  NULL i.e. no local file to use. NOTE if passing a local chain file
  make sure to specify the path to convert from and to the correct build
  like GRCh37 to GRCh38. We can not sense check this for local files.
  The chain file can be submitted as a gz file (as downloaed from
  source) or unzipped.

- rmv_chr:

  Chromosomes to exclude from the formatted summary statistics file. Use
  NULL if no filtering is necessary. Default is
  [`c()`](https://rdrr.io/r/base/c.html) so no chromosomes are removed
  (SNPs on nonstandard chromosomes are always removed). Use
  `c("X", "Y", "MT")` to remove all non-autosomal SNPs.

- verbose:

  Print messages.

## Value

Lifted summary stats in `data.table` or
[GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
format.

## Examples

``` r
sumstats_dt <- MungeSumstats::formatted_example()
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
#> Sorting coordinates with 'data.table'.

sumstats_dt_hg38 <- liftover(sumstats_dt=sumstats_dt, 
                             ref_genome = "hg19",
                             convert_ref_genome="hg38")
#> Performing data liftover from hg19 to hg38.
#> First - Checking CHR format to avoid liftover issues with chain file.
#> Formatting CHR column (if necessary)- Removing chr prefix, converting 23/24 to X/Y and M to MT
#> Converting summary statistics to GenomicRanges.
#> Downloading chain file...
#> Downloading chain file from Ensembl.
#> /tmp/Rtmpt871Bb/GRCh37_to_GRCh38.chain.gz
#> Reordering so first three column headers are SNP, CHR and BP in this order.
#> Reordering so the fourth and fifth columns are A1 and A2.
```
