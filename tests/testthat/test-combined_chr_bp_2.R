test_that("Can correctly separate two CHR:BP columns", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and combine CHR BP columns
  sumstats_dt <- data.table::fread(file)
  #Keep Org to validate values
  sumstats_dt_missing <- data.table::copy(sumstats_dt)
  sumstats_dt_missing[,CHR_BP:=paste0(CHR,":",POS)]
  sumstats_dt_missing[,CHR_BP_2:=paste0(CHR,":",POS)]
  sumstats_dt_missing[,CHR:=NULL]
  sumstats_dt_missing[,POS:=NULL]
  data.table::fwrite(x=sumstats_dt_missing, file=file, sep="\t")
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37")
  res_dt <- data.table::fread(reformatted)
  #Should give same result as separated
  data.table::fwrite(x=sumstats_dt, file=file, sep="\t")
  org <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37")
  org_dt <- data.table::fread(org)
  expect_equal(org_dt,res_dt)
})
