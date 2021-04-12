test_that("Input has correct format", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file)
  reformatted_lines <- readLines(reformatted)
  #Only issue with eduAttainOkbay is the SNP ID name so update before check
  len_eduAttainOkbay <- length(MungeSumstats::eduAttainOkbay)
  expect_equal(reformatted_lines[2:length(reformatted_lines)],
                MungeSumstats::eduAttainOkbay[2:len_eduAttainOkbay])
  headers_eduAttainOkbay <- "SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP"
  expect_equal(reformatted_lines[1],headers_eduAttainOkbay)
})
