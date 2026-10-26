# Infer genome builds

Infers the genome build of summary statistics files (GRCh37 or GRCh38)
from the data. Uses SNP (RSID) & CHR & BP to get genome build.

## Usage

``` r
get_genome_builds(
  sumstats_list,
  header_only = TRUE,
  sampled_snps = 10000,
  names_from_paths = FALSE,
  dbSNP = 155,
  dbSNP_tarball = NULL,
  nThread = 1,
  chr_filt = NULL
)
```

## Arguments

- sumstats_list:

  A named list of paths to summary statistics, or a named list of
  `data.table` objects.

- header_only:

  Instead of reading in the entire `sumstats` file, only read in the
  first N rows where N=`sampled_snps`. This should help speed up cases
  where you have to read in `sumstats` from disk each time.

- sampled_snps:

  Downsample the number of SNPs used when inferring genome build to save
  time.

- names_from_paths:

  Infer the name of each item in `sumstats_list` from its respective
  file path. Only works if `sumstats_list` is a list of paths.

- dbSNP:

  version of dbSNP to be used (144 or 155). Default is 155.

- dbSNP_tarball:

  Pass local versions of dbSNP in tarball format. Default of NULL uses
  the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball` was
  enabled to help with dbSNP versions \>=156, after the decision to no
  longer provide dbSNP releases as bioconductor packages. dbSNP 156
  tarball is available here: http://149.165.171.124/SNPlocs/.

- nThread:

  Number of threads to use for parallel processes.

- chr_filt:

  Internal for testing - filter reference genomes and sumstats to
  specific chromosomes for testing. Pass a list of chroms in format:
  c("1","2"). Default is NULL i.e. no filtering

## Value

ref_genome the genome build of the data

## Details

Iterative version of `get_genome_build`.

## Examples

``` r
# Pass path to Educational Attainment Okbay sumstat file to a temp directory

eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
    package = "MungeSumstats"
)
sumstats_list <- list(ss1 = eduAttainOkbayPth, ss2 = eduAttainOkbayPth)

## Call uses reference genome as default with more than 2GB of memory,
## which is more than what 32-bit Windows can handle so remove certain checks
is_32bit_windows <-
    .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
if (!is_32bit_windows) {
    
    #multiple sumstats can be passed at once to get all their genome builds:
    #ref_genomes <- get_genome_builds(sumstats_list = sumstats_list)
    #just passing first here for speed
    sumstats_list_quick <- list(ss1 = eduAttainOkbayPth)
    ref_genomes <- get_genome_builds(sumstats_list = sumstats_list_quick,
                                     dbSNP=144)
}
#> Inferring genome build of 1 sumstats file(s).
#> Inferring genome build.
#> Reading in only the first 10000 rows of sumstats.
#> Importing tabular file: /home/runner/work/_temp/Library/MungeSumstats/extdata/eduAttainOkbay.txt
#> Checking for empty columns.
#> Standardising column headers.
#> First line of summary statistics file: 
#> MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    
#> Loading SNPlocs data for build 144 on GRCH37.
#> Loading reference genome data.
#> Preprocessing RSIDs.
#> Validating RSIDs of 93 SNPs using BSgenome::snpsById...
#> BSgenome::snpsById done in 10 seconds.
#> Loading SNPlocs data for build 144 on GRCH38.
#> Warning: replacing previous import ‘utils::findMatches’ by ‘S4Vectors::findMatches’ when loading ‘SNPlocs.Hsapiens.dbSNP144.GRCh38’
#> Loading reference genome data.
#> Preprocessing RSIDs.
#> Validating RSIDs of 93 SNPs using BSgenome::snpsById...
#> BSgenome::snpsById done in 29 seconds.
#> Inferred genome build: GRCH37
#> Time difference of 42.72219 secs
#> GRCH37: 1 file(s)
```
