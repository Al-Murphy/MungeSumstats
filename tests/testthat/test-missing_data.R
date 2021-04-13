test_that("Handle missing data", {
  file <- tempfile()
  #Remove data from line 3 to check it is deleted
  eduAttainOkbay_missing <- MungeSumstats::eduAttainOkbay
  eduAttainOkbay_missing[3] <-
    "rs12987662\t2\t100821548\tA\tC\t0.3787\t0.027\t0.003\t"
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
