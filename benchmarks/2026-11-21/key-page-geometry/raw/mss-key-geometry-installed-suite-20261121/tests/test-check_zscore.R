test_that("Check that Z-scores get computed correctly", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        path <- system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats")
        sumstats_dt <- MungeSumstats::read_sumstats(path = path)
        sumstats_dt <- MungeSumstats:::check_zscore(
            sumstats_dt = sumstats_dt, imputation_ind = FALSE,
            compute_z = 'P',force_new_z = TRUE,
            standardise_headers = TRUE,
            mapping_file = sumstatsColHeaders
        )[["sumstats_dt"]]
        sumstats_dt[, newZ := sign(BETA) * sqrt(stats::qchisq(P, 1, 
                                                              lower = FALSE))]
        all_z_equal <- all(sumstats_dt$Z == sumstats_dt$newZ)
        expect_equal(all_z_equal, TRUE)
        #also test BETA/SE
        sumstats_dt <- MungeSumstats:::check_zscore(
          sumstats_dt = sumstats_dt, imputation_ind = FALSE,
          compute_z = 'BETA',force_new_z = TRUE,
          standardise_headers = TRUE,
          mapping_file = sumstatsColHeaders
        )[["sumstats_dt"]]
        sumstats_dt[, newZ := BETA/SE]
        all_z_equal <- all(sumstats_dt$Z == sumstats_dt$newZ)
        expect_equal(all_z_equal, TRUE)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})
