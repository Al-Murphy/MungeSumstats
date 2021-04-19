test_that("Can handle space delimited files", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and write as space delimited
  sumstats_dt <- data.table::fread(file)
  data.table::fwrite(x=sumstats_dt, file=file, sep=" ")
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE)
  res_dt <- data.table::fread(reformatted)
  #check against results of normal run should be the exact same
  file2 <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file2)
  org <- MungeSumstats::format_sumstats(file2,ref_genome="GRCh37",
                                        on_ref_genome = FALSE,
                                        strand_ambig_filter=FALSE,
                                        bi_allelic_filter=FALSE,
                                        allele_flip_check=FALSE)
  org_dt <- data.table::fread(org)
  expect_equal(res_dt,org_dt)
})
