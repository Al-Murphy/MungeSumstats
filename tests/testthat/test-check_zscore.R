test_that("Check that Z-scores get computed correctly", {
    path <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
    sumstats_dt <- MungeSumstats::read_sumstats(path = path)
    sumstats_dt <- check_zscore(sumstats_dt=sumstats_dt,
                                standardise_headers=TRUE)[["sumstats_dt"]]
    sumstats_dt[,newZ:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))]  
    all_z_equal <- all(sumstats_dt$Z==sumstats_dt$newZ)
  expect_equal(all_z_equal, TRUE)
})
