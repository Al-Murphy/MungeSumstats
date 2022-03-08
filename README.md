`MungeSumstats`: Standardise the format of GWAS summary statistics
================
<h5>
<i>Authors</i>: Alan Murphy, Brian Schilder and Nathan Skene
</h5>
<h5>
<i>Updated</i>: Mar-08-2022
</h5>

<!-- Readme.md is generated from Readme.Rmd. Please edit that file -->
<!-- badges: start -->

[![](https://img.shields.io/badge/release%20version-1.2.3-black.svg)](https://www.bioconductor.org/packages/MungeSumstats)
[![](https://img.shields.io/badge/devel%20version-1.3.8-black.svg)](https://github.com/neurogenomics/MungeSumstats)
[![R build
status](https://github.com/neurogenomics/MungeSumstats/workflows/DockerHub/badge.svg)](https://github.com/neurogenomics/MungeSumstats/actions)
[![R build
status](https://github.com/neurogenomics/MungeSumstats/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/neurogenomics/MungeSumstats/actions)
[![](https://img.shields.io/github/last-commit/neurogenomics/MungeSumstats.svg)](https://github.com/neurogenomics/MungeSumstats/commits/master)
[![](https://codecov.io/gh/neurogenomics/MungeSumstats/branch/master/graph/badge.svg)](https://codecov.io/gh/neurogenomics/MungeSumstats)
[![](https://img.shields.io/badge/download-721/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/MungeSumstats)
[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/doi-https://doi.org/10.1093/bioinformatics/btab665-blue.svg)](https://doi.org/https://doi.org/10.1093/bioinformatics/btab665)
<!-- badges: end -->

# Introduction

The `MungeSumstats` package is designed to facilitate the
standardisation of GWAS summary statistics as utilised in our Nature
Genetics paper.<sup>1</sup>

## Overview

The package is designed to handle the lack of standardisation of output
files by the GWAS community. The [MRC IEU Open
GWAS](https://gwas.mrcieu.ac.uk/) team have provided full summary
statistics for &gt;10k GWAS, which are API-accessible via the
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
[Bioconductor](https://bioconductor.org/packages/MungeSumstats)
(≥v3.13). To install `MungeSumstats` on Bioconductor run:

    if (!require("BiocManager")) install.packages("BiocManager")

    BiocManager::install("MungeSumstats")

You can then load the package and data package:

``` r
library(MungeSumstats)
```

Note that for a number of the checks implored by `MungeSumstats` a
reference genome is used. If your GWAS summary statistics file of
interest relates to *GRCh38*, you will need to install
`SNPlocs.Hsapiens.dbSNP144.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`
from Bioconductor as follows:

``` r
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
```

If your GWAS summary statistics file of interest relates to *GRCh37*,
you will need to install `SNPlocs.Hsapiens.dbSNP144.GRCh37` and
`BSgenome.Hsapiens.1000genomes.hs37d5` from Bioconductor as follows:

``` r
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
```

These may take some time to install and are not included in the package
as some users may only need one of *GRCh37*/*GRCh38*. If you are unsure
of the genome build, MungeSumstats can also infer this information from
your data.

# Getting started

See the [Getting started vignette
website](https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html)
for up-to-date instructions on usage.

See the [OpenGWAS vignette
website](https://neurogenomics.github.io/MungeSumstats/articles/OpenGWAS.html)
for information on how to use MungeSumstats to access, standardise and
perform quality control on GWAS Summary Statistics from the MRC IEU
[Open GWAS Project](https://gwas.mrcieu.ac.uk/).

If you have any problems please do file an
[Issue](https://github.com/neurogenomics/MungeSumstats/issues) here on
GitHub.

# Future Enhancements

The `MungeSumstats` package aims to be able to handle the most common
summary statistic file formats including VCF. If your file can not be
formatted by `MungeSumstats` feel free to report the
[Issue](https://github.com/neurogenomics/MungeSumstats/issues) on GitHub
along with your summary statistics file header.

We also encourage people to edit the code to resolve their particular
issues too and are happy to incorporate these through pull requests on
github. If your summary statistic file headers are not recognised by
`MungeSumstats` but correspond to one of

    SNP, BP, CHR, A1, A2, P, Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT, N, N_CAS, N_CON, 
    NSTUDY, INFO or FRQ, 

Feel free to update the `data("sumstatsColHeaders")` following the
approach in the *data.R* file and add your mapping. Then use a [Pull
Request](https://github.com/neurogenomics/MungeSumstats/pulls) on GitHub
and we will incorporate this change into the package.

# Contributors

We would like to acknowledge all those who have contributed to
MungeSumstats:

-   [Alan Murphy](https://github.com/Al-Murphy)
-   [Nathan Skene](https://github.com/NathanSkene)
-   [Brian Schilder](https://github.com/bschilder)
-   [Shea Andrews](https://github.com/sjfandrews)

# References

<div id="refs" class="references csl-bib-body" line-spacing="2">

<div id="ref-Skene2018" class="csl-entry">

<span class="csl-left-margin">1. </span><span
class="csl-right-inline">Nathan G. Skene, T. E. B., Julien Bryois.
Genetic identification of brain cell types underlying schizophrenia.
*Nature Genetics* (2018).
doi:[10.1038/s41588-018-0129-5](https://doi.org/10.1038/s41588-018-0129-5)</span>

</div>

</div>
