test_that("read_header works", {
    
    is_32bit_windows <-
        .Platform$OS.type == "windows" 
    if (!is_32bit_windows) { 
        #### txt ####
        path <- system.file("extdata", "eduAttainOkbay.txt",
                            package = "MungeSumstats")
        header <- read_header(path)
        testthat::expect_equal(dim(header),c(1,9))
        header2 <- read_header(path, n = NULL)
        testthat::expect_equal(dim(header2),c(93,9))
        
        #### vcf #### 
        path <- system.file("extdata","ALSvcf.vcf", package="MungeSumstats") 
        header3 <- read_header(path = path)
        testthat::expect_length(header3, 528)
        header4 <- read_header(path = path, 
                               skip_vcf_metadata = TRUE)
        testthat::expect_equal(dim(header4),c(2,10))
        
        #### vcf.gz ####
        path_gz <- R.utils::gzip(path, remove=FALSE)
        header5 <- read_header(path = path_gz)
        testthat::expect_length(header5, 528)
        
        #### vcf.bgz ####
        dest <- tempfile(fileext = ".vcf.bgz")
        path_bgz <- Rsamtools::bgzip(path, dest = dest)
        header6 <- read_header(path = path_bgz)
        testthat::expect_length(header6, 528)
        header7 <- read_header(path = path_bgz, 
                               skip_vcf_metadata = TRUE)
        testthat::expect_equal(dim(header7),c(2,10))
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})
