test_that("zzz works", {
    
    library(MungeSumstats)
    opts <- options()
    
    testthat::expect_equal(opts$googleAuthR.batch_endpoint,
                           "https://www.googleapis.com/batch")
    testthat::expect_equal(opts$gargle_oauth_cache,
                           "ieugwasr_oauth")
    testthat::expect_equal(opts$googleAuthR.webapp.client_secret,
                           "I7Gqp83Ku4KJxL9zHWYxG_gD")
})
