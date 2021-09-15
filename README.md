Standardise the format of GWAS summary statistics with *MungeSumstats*
================
Alan Murphy, Brian Schilder and Nathan Skene
2021-09-14

<!-- Readme.md is generated from Readme.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/neurogenomics/MungeSumstats/workflows/R-full/badge.svg)](https://github.com/neurogenomics/MungeSumstats/actions)
[![Codecov test
coverage](https://codecov.io/gh/neurogenomics/MungeSumstats/branch/master/graph/badge.svg)](https://codecov.io/gh/neurogenomics/MungeSumstats?branch=master)
<!-- badges: end -->

# Introduction

The *MungeSumstats* package is designed to facilitate the
standardisation of GWAS summary statistics as utilised in our Nature
Genetics paper.<sup>1</sup> If you use MungeSumstats, please cite our
preprint [Murphy and Skene. MungeSumstats: A Bioconductor package for
the standardisation and quality control of many GWAS summary
statistics](https://www.biorxiv.org/content/10.1101/2021.06.21.449239v1).

# Overview

The package is designed to handle the lack of standardisation of output
files by the GWAS community. There is a group who have now manually
standardised many GWAS: [R interface to the IEU GWAS database API •
ieugwasr](https://mrcieu.github.io/ieugwasr/) and
[gwasvcf](https://github.com/MRCIEU/gwasvcf) but because a lot of GWAS
remain closed access, these repositories are not all encompassing.

*MungeSumstats* provides a framework to standardise the format for any
GWAS summary statistics, including those in VCF format, enabling
downstream integration and analysis. The package works by addressing the
most common discrepancies across summary statistic files.
*MungeSumstats* also offers a range of adjustable, Quality Control (QC).

# Installing MungeSumstats

*MungeSumstats* is avaiable on
[Bioconductor](https://bioconductor.org/packages/MungeSumstats) (v3.13).
To install *MungeSumstats* on Bioconductor run:

    if (!require("BiocManager"))
        install.packages("BiocManager")

    BiocManager::install("MungeSumstats")

You can then load the package and data package:

``` r
library(MungeSumstats)
```

Note that for a number of the checks implored by *MungeSumstats* a
reference genome is used. If your GWAS summary statistics file of
interest relates to *GRCh38*, you will need to install
`SNPlocs.Hsapiens.dbSNP144.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`
from Bioconductor as follows:

    BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

If your GWAS summary statistics file of interest relates to *GRCh37*,
you will need to install `SNPlocs.Hsapiens.dbSNP144.GRCh37` and
`BSgenome.Hsapiens.1000genomes.hs37d5` from Bioconductor as follows:

    BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

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

If you have any problems please do file an issue here on github.

# Citation

If you use the MungeSumstats package then please cite

[Murphy and Skene. MungeSumstats: A Bioconductor package for the
standardisation and quality control of many GWAS summary
statistics](https://www.biorxiv.org/content/10.1101/2021.06.21.449239v1).

# Future Enhancements

The *MungeSumstats* package aims to be able to handle the most common
summary statistic file formats including VCF. If your file can not be
formatted by *MungeSumstats* feel free to report the bug on github:
<https://github.com/neurogenomics/MungeSumstats> along with your summary
statistic file header.

We also encourage people to edit the code to resolve their particular
issues too and are happy to incorporate these through pull requests on
github. If your summary statistic file headers are not recognised by
*MungeSumstats* but correspond to one of

    SNP, BP, CHR, A1, A2, P, Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT, N, N_CAS, N_CON, 
    NSTUDY, INFO or FRQ, 

feel free to update the `data("sumstatsColHeaders")` following the
approach in the data.R file and add your mapping. Then use a pull
request on github and we will incorporate this change into the package.

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
