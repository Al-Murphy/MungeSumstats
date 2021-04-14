test_that("validate parameters", {
  
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  
  #Check that submitting incorrect parameters gives errors
  #Run MungeSumstats code
  #passing incorrect path
  not_file <- tempfile()
  error_return <-
    tryCatch( MungeSumstats::format_sumstats(not_file,ref_genome="GRCh37"),
              error = function(e) e,
              warning = function(w) w
    )
  #incorrect ref genome
  error_return2 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="fake"),
              error = function(e) e,
              warning = function(w) w
    )
  #incorrect convert_small_p value
  error_return3 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                              convert_small_p = "fake"),
              error = function(e) e,
              warning = function(w) w
    )
  #incorrect convert_n_int value
  error_return4 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             convert_n_int = "fake"),
              error = function(e) e,
              warning = function(w) w
    )
  expect_equal(is(error_return, "error"), TRUE)
  expect_equal(is(error_return2, "error"), TRUE)
  expect_equal(is(error_return3, "error"), TRUE)
  expect_equal(is(error_return4, "error"), TRUE)
})
