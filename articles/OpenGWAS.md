# Import GWAS summary statistics from Open GWAS

``` r

library(MungeSumstats)
```

MungeSumstats now offers high throughput query and import functionality
to data from the MRC IEU [Open GWAS
Project](https://gwas.mrcieu.ac.uk/).

This is made possible by the use the [IEU OpwnGWAS R
package](https://github.com/MRCIEU/ieugwasr): `ieugwasr`.

Before you can use this functionality however, please complete the
following steps:

## Authenticate Access IEU OpenGWAS API

To authenticate, you need to generate a token from the OpenGWAS website.
The token behaves like a password, and it will be used to authorise the
requests you make to the OpenGWAS API. Here are the steps to generate
the token and then have `ieugwasr` automatically use it for your
queries:

1.  Login to <https://api.opengwas.io/profile/>
2.  Generate a new token
3.  Add `OPENGWAS_JWT=<token>` to your .Renviron file, thi can be edited
    in R by running
    [`usethis::edit_r_environ()`](https://usethis.r-lib.org/reference/edit.html)
4.  Restart your R session
5.  To check that your token is being recognised, run
    [`ieugwasr::get_opengwas_jwt()`](https://mrcieu.github.io/ieugwasr/reference/get_opengwas_jwt.html).
    If it returns a long random string then you are authenticated.
6.  To check that your token is working, run
    [`ieugwasr::user()`](https://mrcieu.github.io/ieugwasr/reference/user.html).
    It will make a request to the API for your user information using
    your token. It should return a list with your user information. If
    it returns an error, then your token is not working.
7.  Make sure you have submitted user information to increasse you API
    limit at <https://api.opengwas.io/profile/>.

## Find GWAS datasets

We can search by terms and with other filters like sample size:

``` r

#### Search for datasets ####
metagwas <- MungeSumstats::find_sumstats(traits = c("parkinson","alzheimer"), 
                                         min_sample_size = 1000)
head(metagwas,3)
ids <- (dplyr::arrange(metagwas, nsnp))$id  
```

    ##          id               trait group_name year    author
    ## 1 ieu-a-298 Alzheimer's disease     public 2013   Lambert
    ## 2   ieu-b-2 Alzheimer's disease     public 2019 Kunkle BW
    ## 3 ieu-a-297 Alzheimer's disease     public 2013   Lambert
    ##                                                                                                                                                                                                                                                                                                                    consortium
    ## 1                                                                                                                                                                                                                                                                                                                        IGAP
    ## 2 Alzheimer Disease Genetics Consortium (ADGC), European Alzheimer's Disease Initiative (EADI), Cohorts for Heart and Aging Research in Genomic Epidemiology Consortium (CHARGE), Genetic and Environmental Risk in AD/Defining Genetic, Polygenic and Environmental Risk for Alzheimer's Disease Consortium (GERAD/PERADES),
    ## 3                                                                                                                                                                                                                                                                                                                        IGAP
    ##                 sex population     unit     nsnp sample_size       build
    ## 1 Males and Females   European log odds    11633       74046 HG19/GRCh37
    ## 2 Males and Females   European       NA 10528610       63926 HG19/GRCh37
    ## 3 Males and Females   European log odds  7055882       54162 HG19/GRCh37
    ##   category                subcategory ontology mr priority     pmid sd
    ## 1  Disease Psychiatric / neurological       NA  1        1 24162737 NA
    ## 2   Binary Psychiatric / neurological       NA  1        0 30820047 NA
    ## 3  Disease Psychiatric / neurological       NA  1        2 24162737 NA
    ##                                                                      note ncase
    ## 1 Exposure only; Effect allele frequencies are missing; forward(+) strand 25580
    ## 2                                                                      NA 21982
    ## 3                Effect allele frequencies are missing; forward(+) strand 17008
    ##   ncontrol     N
    ## 1    48466 74046
    ## 2    41944 63926
    ## 3    37154 54162

You can also search by ID:

``` r

### By ID and sample size
metagwas <- find_sumstats(
  ids = c("ieu-b-4760", "prot-a-1725", "prot-a-664"),
  min_sample_size = 5000
)
```

## Import full results

You can supply
[`import_sumstats()`](https://al-murphy.github.io/MungeSumstats/reference/import_sumstats.md)
with a list of as many OpenGWAS IDs as you want, but we’ll just give one
to save time.

``` r

datasets <- MungeSumstats::import_sumstats(ids = "ieu-a-298",
                                           ref_genome = "GRCH37")
```

### Summarise results

By default, `import_sumstats` results a named list where the names are
the Open GWAS dataset IDs and the items are the respective paths to the
formatted summary statistics.

``` r

print(datasets)
```

    ## $`ieu-a-298`
    ## [1] "/tmp/RtmpxsovP7/ieu-a-298.tsv.gz"

You can easily turn this into a data.frame as well.

``` r

results_df <- data.frame(id=names(datasets), 
                         path=unlist(datasets))
print(results_df)
```

    ##                  id                             path
    ## ieu-a-298 ieu-a-298 /tmp/RtmpxsovP7/ieu-a-298.tsv.gz

## Import full results (parallel)

*Optional*: Speed up with multi-threaded download via
[axel](https://github.com/axel-download-accelerator/axel).

``` r

datasets <- MungeSumstats::import_sumstats(ids = ids, 
                                           vcf_download = TRUE, 
                                           download_method = "axel", 
                                           nThread = max(2,future::availableCores()-2))
```

## Further functionality

See the [Getting started
vignette](https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html)
for more information on how to use MungeSumstats and its functionality.

## Session Info

``` r

utils::sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] MungeSumstats_1.21.2 BiocStyle_2.40.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1            dplyr_1.2.1                
    ##  [3] blob_1.3.0                  R.utils_2.13.0             
    ##  [5] Biostrings_2.80.2           bitops_1.1-0               
    ##  [7] fastmap_1.2.0               RCurl_1.98-1.20            
    ##  [9] VariantAnnotation_1.58.0    GenomicAlignments_1.48.0   
    ## [11] XML_3.99-0.25               digest_0.6.39              
    ## [13] lifecycle_1.0.5             KEGGREST_1.52.2            
    ## [15] RSQLite_3.53.3              magrittr_2.0.5             
    ## [17] compiler_4.6.1              rlang_1.3.0                
    ## [19] sass_0.4.10                 tools_4.6.1                
    ## [21] yaml_2.3.12                 data.table_1.18.6.1        
    ## [23] rtracklayer_1.72.0          knitr_1.52                 
    ## [25] S4Arrays_1.12.1             bit_4.6.0                  
    ## [27] curl_8.1.0                  DelayedArray_0.38.2        
    ## [29] ieugwasr_1.2.0              abind_1.4-8                
    ## [31] BiocParallel_1.46.0         BiocGenerics_0.58.1        
    ## [33] desc_1.4.3                  R.oo_1.27.1                
    ## [35] grid_4.6.1                  stats4_4.6.1               
    ## [37] SummarizedExperiment_1.42.0 cli_3.6.6                  
    ## [39] rmarkdown_2.32              crayon_1.5.3               
    ## [41] ragg_1.5.2                  generics_0.1.4             
    ## [43] otel_0.2.0                  httr_1.4.9                 
    ## [45] rjson_0.2.23                DBI_1.3.0                  
    ## [47] cachem_1.1.0                stringr_1.6.0              
    ## [49] parallel_4.6.1              AnnotationDbi_1.74.0       
    ## [51] BiocManager_1.30.27         XVector_0.52.0             
    ## [53] restfulr_0.0.17             matrixStats_1.5.0          
    ## [55] vctrs_0.7.3                 Matrix_1.7-5               
    ## [57] jsonlite_2.0.0              bookdown_0.48              
    ## [59] IRanges_2.46.0              S4Vectors_0.50.3           
    ## [61] bit64_4.8.6                 systemfonts_1.3.2          
    ## [63] GenomicFeatures_1.64.0      jquerylib_0.1.4            
    ## [65] glue_1.8.1                  pkgdown_2.2.1              
    ## [67] codetools_0.2-20            stringi_1.8.9              
    ## [69] GenomeInfoDb_1.48.0         BiocIO_1.22.0              
    ## [71] GenomicRanges_1.64.0        UCSC.utils_1.8.0           
    ## [73] tibble_3.3.1                pillar_1.11.1              
    ## [75] htmltools_0.5.9             Seqinfo_1.2.0              
    ## [77] BSgenome_1.80.0             R6_2.6.1                   
    ## [79] textshaping_1.0.5           evaluate_1.0.5             
    ## [81] lattice_0.22-9              Biobase_2.72.0             
    ## [83] R.methodsS3_1.8.2           png_0.1-9                  
    ## [85] Rsamtools_2.28.0            cigarillo_1.2.1            
    ## [87] memoise_2.0.1               bslib_0.12.0               
    ## [89] SparseArray_1.12.3          xfun_0.62                  
    ## [91] fs_2.1.0                    MatrixGenerics_1.24.0      
    ## [93] pkgconfig_2.0.3
