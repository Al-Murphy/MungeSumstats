test_that("zzz works", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
    ##.Platform$r_arch == "i386"
    #only run test if user has internet access
    #capture internet outage or server issues
    if(try(is.character(getURL("www.google.com")))==TRUE &&
       !is_32bit_windows){
        library(MungeSumstats)
        opts <- options()
        
        #testthat::expect_equal(opts$googleAuthR.batch_endpoint,
        #                       "https://www.googleapis.com/batch")
        testthat::expect_equal(opts$gargle_oauth_cache,
                               "ieugwasr_oauth")
        testthat::expect_equal(opts$googleAuthR.webapp.client_secret,
                               "I7Gqp83Ku4KJxL9zHWYxG_gD")
    }
    else{
        #expect_equal(TRUE, TRUE)
        expect_equal(TRUE, TRUE)
        expect_equal(TRUE, TRUE)
    }
})
