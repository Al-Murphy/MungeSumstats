test_that("Handle n when it isn't an integer", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and make N
  sumstats_dt <- data.table::fread(file)
  #Add N column and make it not an integer
  sumstats_dt[,N:=10*runif(nrow(sumstats_dt))]
  sumstats_dt[,N_fixed:=round(N,0)]
  data.table::fwrite(x=sumstats_dt, file=file, sep="\t")
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE)
  #In results if N = N_fixed it worked
  res_dt <- data.table::fread(reformatted)
  expect_equal(res_dt$N,res_dt$N_FIXED)
})
