Standardise the format of summary statistics from GWAS with
*MungeSumstats*
================
Alan Murphy and Nathan Skene
2021-04-28

<!-- Readme.md is generated from Readme.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/neurogenomics/MungeSumstats/workflows/R-full/badge.svg)](https://github.com/neurogenomics/MungeSumstats/actions)
[![Codecov test
coverage](https://codecov.io/gh/neurogenomics/MungeSumstats/branch/master/graph/badge.svg)](https://codecov.io/gh/neurogenomics/MungeSumstats?branch=master)
<!-- badges: end -->

# Introduction

The *MungeSumstats* package is designed to facilitate the
standardisation of GWAS summary statistics as utilised in our Nature
Genetics paper.<sup>1</sup>

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
*MungeSumstats* also offers a range of adjustable, Quality Control (QC)
steps.

# Installing MungeSumstats

To install *MungeSumstats* on Bioconductor run:

    if (!require("BiocManager"))
        install.packages("BiocManager")

    BiocManager::install("MungeSumstats")

You can then load the package and data package:

``` r
library(MungeSumstats)
```

# Getting started

See the vignette for use cases of *MungeSumstats*:

    browseVignettes("MungeSumstats")

If you have any problems please do file an issue here on github.

# Citation

If you use the MungeSumstats package as well then please cite

[Skene, et al. Genetic identification of brain cell types underlying
schizophrenia. Nature Genetics,
2018.](https://www.nature.com/articles/s41588-018-0129-5)

# Future Enhancements

The *MungeSumstats* package aims to be able to handle the most common
summary statistic file formats including VCF. If your file can not be
formatted by *MungeSumstats* feel free to report the bug on github:
<https://github.com/neurogenomics/MungeSumstats> along with your summary
statistic file header.

We also encourage people to edit the code to resolve their particular
issues too and are happy to incorporate these through pull requests on
github. If your summary statistic file headers are not recognised by
*MungeSumstats* but correspond to one of SNP, BP, CHR, A1, A2, P, Z, OR,
BETA, LOG\_ODDS, SIGNED\_SUMSTAT, N, N\_CAS, N\_CON, NSTUDY, INFO or
FRQ, feel free to update the `MungeSumstats::sumstatsColHeaders`
following the approach in the data.R file and add your mapping. Then use
a pull request on github and we will incorporate this change into the
package.

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
