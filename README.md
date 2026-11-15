`MungeSumstats`: Standardise the format of GWAS summary statistics
================
<h5>  
<i>Authors</i>: Alan Murphy, Brian Schilder and Nathan Skene  
</h5>
<h5>  
<i>Updated</i>: Nov-06-2025  
</h5>

<!-- Readme.md is generated from Readme.Rmd. Please edit that file -->

<!-- badges: start -->

[![](https://img.shields.io/badge/release%20version-1.16.0-black.svg)](https://www.bioconductor.org/packages/MungeSumstats)
[![](https://img.shields.io/badge/devel%20version-1.13.4-black.svg)](https://github.com/neurogenomics/MungeSumstats)
[![](https://img.shields.io/github/last-commit/neurogenomics/MungeSumstats.svg)](https://github.com/neurogenomics/MungeSumstats/commits/master)
[![](https://img.shields.io/badge/download-31127/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/MungeSumstats)
[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/doi-https://doi.org/10.1093/bioinformatics/btab665-blue.svg)](https://doi.org/https://doi.org/10.1093/bioinformatics/btab665)
<!-- badges: end -->
<!--[![R build status](https://github.com/neurogenomics/MungeSumstats/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/MungeSumstats/actions)-->
<!--`[![](https://codecov.io/gh/neurogenomics/MungeSumstats/branch/master/graph/badge.svg)](https://app.codecov.io/gh/neurogenomics/MungeSumstats) -->

# Introduction

The `MungeSumstats` package is designed to facilitate the
standardisation of GWAS summary statistics.

## Overview

The package is designed to handle the lack of standardisation of output
files by the GWAS community. The [MRC IEU Open
GWAS](https://gwas.mrcieu.ac.uk/) team have provided full summary
statistics for \>10k GWAS, which are API-accessible via the
[`ieugwasr`](https://mrcieu.github.io/ieugwasr/) and
[`gwasvcf`](https://github.com/MRCIEU/gwasvcf) packages. But these GWAS
are only standardised in the sense that they are VCF format, and can be
fully standardised with `MungeSumstats`.

`MungeSumstats` provides a framework to standardise the format for any
GWAS summary statistics, including those in VCF format, enabling
downstream integration and analysis. It addresses the most common
discrepancies across summary statistic files, and offers a range of
adjustable Quality Control (QC) steps.

## Citation

If you use `MungeSumstats`, please cite the original authors of the GWAS
as well as:

> Alan E Murphy, Brian M Schilder, Nathan G Skene (2021) MungeSumstats:
> A Bioconductor package for the standardisation and quality control of
> many GWAS summary statistics. *Bioinformatics*, btab665,
> <https://doi.org/10.1093/bioinformatics/btab665>

# Installing `MungeSumstats`

`MungeSumstats` is available on
[Bioconductor](https://bioconductor.org/packages/MungeSumstats). To
install `MungeSumstats` on Bioconductor run:

``` r
if (!require("BiocManager")) install.packages("BiocManager")

BiocManager::install("MungeSumstats")
```

You can then load the package and data package:

``` r
library(MungeSumstats)
```

Note that there is also a [docker image for
MungeSumstats](https://hub.docker.com/r/neurogenomicslab/mungesumstats).

Note that for a number of the checks implored by `MungeSumstats` a
reference genome is used. If your GWAS summary statistics file of
interest relates to *GRCh38*, you will need to install
`SNPlocs.Hsapiens.dbSNP155.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`
from Bioconductor as follows:

``` r
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
```

If your GWAS summary statistics file of interest relates to *GRCh37*,
you will need to install `SNPlocs.Hsapiens.dbSNP155.GRCh37` and
`BSgenome.Hsapiens.1000genomes.hs37d5` from Bioconductor as follows:

``` r
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
```

These may take some time to install and are not included in the package
as some users may only need one of *GRCh37*/*GRCh38*. If you are unsure
of the genome build, MungeSumstats can also infer this information from
your data.

# Getting started

See the [Getting started vignette
website](https://al-murphy.github.io/MungeSumstats/articles/MungeSumstats.html)
for up-to-date instructions on usage.

See the [OpenGWAS vignette
website](https://al-murphy.github.io/MungeSumstats/articles/OpenGWAS.html)
for information on how to use MungeSumstats to access, standardise and
perform quality control on GWAS Summary Statistics from the MRC IEU
[Open GWAS Project](https://gwas.mrcieu.ac.uk/).

**NOTE** to authenticate, you need to generate a token from the OpenGWAS
website. The token behaves like a password, and it will be used to
authorise the requests you make to the OpenGWAS API. Here are the steps
to generate the token and then have `ieugwasr` automatically use it for
your queries:

1.  Login to <https://api.opengwas.io/profile/>
2.  Generate a new token
3.  Add `OPENGWAS_JWT=<token>` to your .Renviron file, thi can be edited
    in R by running `usethis::edit_r_environ()`
4.  Restart your R session
5.  To check that your token is being recognised, run
    `ieugwasr::get_opengwas_jwt()`. If it returns a long random string
    then you are authenticated.
6.  To check that your token is working, run `ieugwasr::user()`. It will
    make a request to the API for your user information using your
    token. It should return a list with your user information. If it
    returns an error, then your token is not working.
7.  Make sure you have submitted use

Please read carefully through the [FAQ
website](https://github.com/Al-Murphy/MungeSumstats/wiki/FAQ) for an
queries about running MungeSumstats. If you have any outside of this
problems please do file an
[Issue](https://github.com/al-murphy/MungeSumstats/issues) here on
GitHub.

# Future Enhancements

The `MungeSumstats` package aims to be able to handle the most common
summary statistic file formats including VCF. If your file can not be
formatted by `MungeSumstats` feel free to report the
[Issue](https://github.com/al-murphy/MungeSumstats/issues) on GitHub
along with your summary statistics file header.

We also encourage people to edit the code to resolve their particular
issues too and are happy to incorporate these through pull requests on
github. If your summary statistic file headers are not recognised by
`MungeSumstats` but correspond to one of

    SNP, BP, CHR, A1, A2, P, Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT, N, N_CAS, N_CON, 
    NSTUDY, INFO or FRQ, 

Feel free to update the `data("sumstatsColHeaders")` following the
approach in the *data.R* file and add your mapping. Then use a [Pull
Request](https://github.com/al-murphy/MungeSumstats/pulls) on GitHub and
we will incorporate this change into the package.

# Contributors

We would like to acknowledge all those who have contributed to
`MungeSumstats` development:

- [Alan Murphy](https://github.com/Al-Murphy)
- [Nathan Skene](https://github.com/NathanSkene)
- [Brian Schilder](https://github.com/bschilder)
- [Shea Andrews](https://github.com/sjfandrews)
- [Jonathan Griffiths](https://github.com/jonathangriffiths)
- [Kitty Murphy](https://github.com/KittyMurphy)
- [Mykhaylo Malakhov](https://github.com/MykMal)
- [Alasdair Warwick](https://github.com/rmgpanw)
- [Ao Lu](https://github.com/leoarrow1)
- [Sufyan Sualeman](https://github.com/sufyansuleman)
- [Toby Clark](https://github.com/toby-clark4)

### Preservation-oriented portable interface

The separate `prepare_reference()`, `harmonise_sumstats()` and `harmonise_file()`
APIs preserve every source row, identifier and original measurement token. The
existing `format_sumstats()` filtering behavior is unchanged. A C++17 core trims
and left-aligns sequence-resolved SNP/indel contrasts against a declared plain
FASTA, with disk-backed one-to-many coordinate/ID/merged-alias lookup. Runtime
normalization requires no external bcftools or Rhtslib.

Reference TSVs have exactly `chr`, `pos`, `ref`, `alt`, `id`, `aliases` columns in
that order. Positions are 1-based, ALT can contain comma-separated alleles, and
aliases are semicolon separated. Contig names must match the FASTA exactly.

```r
reference <- prepare_reference(
  "GRCh38.fa", "reference-variants.tsv", "references/dbsnp-release-GRCh38",
  build = "GRCh38", version = "declared-reference-release")

# A1 = other allele, A2 = effect allele; source remains exactly as supplied.
result <- harmonise_sumstats(sumstats, reference, frequency_type = "EAF")
result$normalized  # nullable, row-aligned native coordinates/effects/status
result$source      # original values, including raw P-value tokens

# Bounded file-to-file path, including gzip input; output path must be new.
handle <- open_reference(reference)  # Verify once, reuse across files/batches.
manifest <- harmonise_file("study.tsv.gz", handle, "editions/study-v1",
                           batch_size = 10000L)
```

`match`, `swap`, `complement` and `swapped_complement` are explicit. Palindromic
or multiple compatible mappings remain ambiguous. ALT–ALT contrasts, missing
rsIDs, unsupported/symbolic alleles and unknown optional fields are preserved;
normalization may remain null. No first ALT is chosen. Swapped EAF is never
computed as `1-p` when a third allele is known, including across split reference
rows. Frequency type defaults to unknown; only an explicit `frequency_type="EAF"`
permits normalized frequencies. MAF/unknown frequencies are not treated as EAF. Liftover is not attempted,
and GRCh38 projection fields remain null.

The file interface supports unquoted TSV/CSV with one record per line and a final
newline. It stages the original artifact, raw source lines, normalized TSV with
`source_row` and `normalized_build`, audit counts and manifest before publishing a new edition directory.
Portable quoted TSV receipts (`metadata.tsv`, `columns.tsv`,
`reference_checksums.tsv`, `checksums.tsv`) support non-R consumers; the checksum
list covers the other artifacts and receipts, excluding itself and `manifest.rds`.
Malformed/truncated reads fail without publication. No parser can identify a
syntactically valid partial acquisition without an external receipt.

The reference sorter uses bounded chunks and merge fan-in of 32; lookups binary
search disk offsets and cap candidates rather than loading reference variants in
RAM. The table interface keeps its source and result in R; use `harmonise_file()`
for a bounded-batch corpus path. FASTA must be plain and fixed-width (except final
sequence lines); unanchored `I/D/-` alleles and compressed reference FASTA are not
interpreted. Source FASTA remains an external immutable file identified by path
and MD5. Hash verification reads the reference sequentially once when opening a live
handle; subsequent calls check file sizes and modification times. Temporary disk requirements include the input snapshot and external sort
runs. No million-row/reference-scale performance, 8-GiB RSS target, biological
recovery rate or Windows execution has yet been validated. A three-platform CI
workflow is supplied to run the standalone fixtures.

Standalone verification does not require the package's Bioconductor imports:

```sh
# Rcpp and testthat must be installed in the active R library.
Rscript tests/standalone/preservation.R
Rscript tests/standalone/registered.R
./tests/native/run.sh
```

Optional independent normalization oracle: `python3 tests/native/oracle.py`
(requires bcftools for testing only).

For caller-supplied benchmark assets, run the file path with an explicit expected
row count and retain the output edition plus metrics:

```sh
Rscript scripts/benchmark-preservation.R \
  --input study.tsv.gz --reference references/dbsnp-release-GRCh38 \
  --output benchmarks/study-batch10000 --metrics benchmarks/study-batch10000.tsv \
  --rows 1000000 --batch-size 10000 --frequency-type EAF --effect-type beta \
  --columns chr=CHR,pos=BP,id=SNP,other=A1,effect_allele=A2,effect=BETA,eaf=FRQ \
  --cache-state unspecified-uncontrolled
```

The harness separately times reference opening/verification and file processing,
checks source/output counts against `--rows`, and records provenance and hashes.
`--source-root /path/to/MungeSumstats` supports standalone Rcpp compilation when
full package dependencies are unavailable; compilation is outside the timed
sections. Cache state is caller declared and never flushed or verified by this
script. R allocations are not RSS: collect process memory externally with Linux
`/usr/bin/time -v`, macOS `/usr/bin/time -l`, or a Windows process monitor, retaining
the raw output. External process measurements with `--source-root` can include
compilation; use an installed package for representative runtime measurements.
The tiny `Rscript tests/standalone/benchmark.R` smoke test makes no performance
inference. Reference-scale 1M/10M/40M runs require separately supplied assets.
