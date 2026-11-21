test_that("check_range_p_val works", {
    
    sumstats_dt <- MungeSumstats:::formatted_example()
    range1 <- 1:3
    range2 <- 6:10
    sumstats_dt$P[range1] <- 5
    sumstats_dt$P[range2] <- -5
    range_n <- length(range1) + length(range2)
    sumstats <- check_range_p_val(sumstats_dt = sumstats_dt,
                                  convert_large_p = TRUE,
                                  convert_neg_p = TRUE,
                                  imputation_ind = TRUE)
    testthat::expect_true(
        "convert_large_p_1" %in% colnames(sumstats$sumstats_dt))
    testthat::expect_equal(
        sum(sumstats$sumstats_dt$convert_large_p_1,na.rm = TRUE),
        length(range1))
    
    testthat::expect_true(
        "convert_neg_p_0" %in% colnames(sumstats$sumstats_dt))
    testthat::expect_equal(
        sum(sumstats$sumstats_dt$convert_neg_p_0,na.rm = TRUE),
        length(range2))
})
