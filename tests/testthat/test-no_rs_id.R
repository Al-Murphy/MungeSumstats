test_that("Handle rs ID in row", {
  file <- tempfile()
  #Remove data from line 3 to check it is deleted
  eduAttainOkbay_missing <- MungeSumstats::eduAttainOkbay
  eduAttainOkbay_missing[3] <-
    "9320913\t6\t98584733\tA\tC\t0.5019\t0.024\t0.003\t2.457e-19"
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(eduAttainOkbay_missing,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37")
  reformatted_lines <- readLines(reformatted)
  #Should equal org apart from this one line
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  org <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37")
  org_lines <- readLines(org)
  expect_equal(reformatted_lines,org_lines[-3])
})
