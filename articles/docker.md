# Docker/Singularity Containers

## Installation

MungeSumstats is now available via
[ghcr.io](https://ghcr.io/ghcr.io/neurogenomics/MungeSumstats) as a
containerised environment with Rstudio and all necessary dependencies
pre-installed.

### Method 1: via Docker

First, [install Docker](https://docs.docker.com/get-docker/) if you have
not already.

Create an image of the [Docker](https://www.docker.com/) container in
command line:

    docker pull ghcr.io/neurogenomics/MungeSumstats

Once the image has been created, you can launch it with:

    docker run \
      -d \
      -e ROOT=true \
      -e PASSWORD="<your_password>" \
      -v ~/Desktop:/Desktop \
      -v /Volumes:/Volumes \
      -p 8900:8787 \
      ghcr.io/neurogenomics/MungeSumstats

#### NOTES

- Make sure to replace `<your_password>` above with whatever you want
  your password to be.  
- Change the paths supplied to the `-v` flags for your particular use
  case.
- The `-d` ensures the container will run in “detached” mode, which
  means it will persist even after you’ve closed your command line
  session.  
- The username will be *“rstudio”* by default.  
- Optionally, you can also install the [Docker
  Desktop](https://www.docker.com/products/docker-desktop/) to easily
  manage your containers.

### Method 2: via Singularity

If you are using a system that does not allow Docker (as is the case for
many institutional computing clusters), you can instead [install Docker
images via
Singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html).

    singularity pull docker://ghcr.io/neurogenomics/MungeSumstats

For troubleshooting, see the [Singularity
documentation](https://docs.sylabs.io/guides/latest/user-guide/singularity_and_docker.html#github-container-registry).

## Usage

Finally, launch the containerised Rstudio by entering the following URL
in any web browser: *<http://localhost:8900/>*

Login using the credentials set during the Installation steps.

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

  
