test_that("Handle duplicate rows based on rs ID", {
  file <- tempfile()
  #Duplicate rows and check output is the same
  eduAttainOkbay_missing <- MungeSumstats::eduAttainOkbay
  len <- length(eduAttainOkbay_missing)
  eduAttainOkbay_missing <- c(eduAttainOkbay_missing,
                                eduAttainOkbay_missing[2:len])
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(eduAttainOkbay_missing,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37")
  reformatted_lines <- readLines(reformatted)
  #Should equal org apart from this one line
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  org <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37")
  org_lines <- readLines(org)
  #check equal regardless of order
  expect_equal(setequal(reformatted_lines,org_lines),TRUE)
})
