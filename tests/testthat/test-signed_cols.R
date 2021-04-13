test_that("will fail without a signed column", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and combine CHR BP columns
  sumstats_dt <- data.table::fread(file)
  #Remove signed column - BETA
  sumstats_dt[,Beta:=NULL]
  data.table::fwrite(x=sumstats_dt, file=file, sep="\t")
  #Run MungeSumstats code
  error_return <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37"),
              error = function(e) e,
              warning = function(w) w
    )
  expect_equal(is(error_return, "error"), TRUE)
})
