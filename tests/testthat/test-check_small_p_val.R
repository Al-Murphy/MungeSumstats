test_that("check_small_p_val works", {
  
    sumstats_dt <- MungeSumstats:::formatted_example()
    range1 <- 1:3
    range2 <- 6:10
    range_n <- length(range1) + length(range2)
    sumstats_dt$P[range1] <- 5e-324
    sumstats_dt$P[range2] <- "5e-324"
    sumstats <- check_small_p_val(sumstats_dt = sumstats_dt,
                                  convert_small_p = TRUE,
                                  imputation_ind = TRUE)
    testthat::expect_true(
        "convert_small_p_0" %in% colnames(sumstats$sumstats_dt))
    testthat::expect_equal(
        sum(sumstats$sumstats_dt$convert_small_p_0,na.rm = TRUE),
        range_n)
})
