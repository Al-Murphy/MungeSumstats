test_that("Input has correct headers", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and correct column headers
  sumstats_dt <- data.table::fread(file)
  names(sumstats_dt) <- c("SNP","CHR","BP","A1","A2","FRQ","BETA","SE","P")
  data.table::fwrite(x=sumstats_dt, file=file, sep="\t")
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE)
  reformatted_dt <- data.table::fread(reformatted)
  expect_equal(all.equal(reformatted_dt,sumstats_dt,ignore.row.order=TRUE),
                TRUE)
})
