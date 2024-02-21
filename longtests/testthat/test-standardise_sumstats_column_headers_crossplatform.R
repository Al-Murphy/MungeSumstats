test_that("standardise_sumstats_column_headers_crossplatform works", {
    
    sumstats_dt <- data.table::fread(
        system.file("extdata", "eduAttainOkbay.txt",
                    package = "MungeSumstats"))
    sumstats_dt$Support <- 1
    #### Internal version ####
    sumstats_dt2 <- 
        MungeSumstats:::standardise_sumstats_column_headers_crossplatform(
            ## Important! Must copy input data to avoid it also being modified 
            ## by the function.
            sumstats_dt = data.table::copy(sumstats_dt), 
            uppercase_unmapped = FALSE,
            return_list = FALSE)
    testthat::expect_true("Support" %in% colnames(sumstats_dt))
    testthat::expect_true("Pval" %in% colnames(sumstats_dt))
    testthat::expect_true("Support" %in% colnames(sumstats_dt2))
    testthat::expect_true("P" %in% colnames(sumstats_dt2))
    
    
    #### Exported version ####
    sumstats_dt3 <- 
        MungeSumstats::standardise_header(
            ## Important! Must copy input data to avoid it also being modified 
            ## by the function.
            sumstats_dt = data.table::copy(sumstats_dt),  
            uppercase_unmapped = FALSE,
            return_list = FALSE)
    testthat::expect_true("Support" %in% colnames(sumstats_dt))
    testthat::expect_true("Pval" %in% colnames(sumstats_dt))
    testthat::expect_true("Support" %in% colnames(sumstats_dt3))
    testthat::expect_true("P" %in% colnames(sumstats_dt3))
})
