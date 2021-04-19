test_that("Imputation of CHR and BP correctly", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and drop CHR BP columns
  sumstats_dt <- data.table::fread(file)
  #Keep Org to validate values
  sumstats_dt_missing <- data.table::copy(sumstats_dt)
  sumstats_dt_missing[,CHR:=NULL]
  sumstats_dt_missing[,POS:=NULL]
  data.table::fwrite(x=sumstats_dt_missing, file=file, sep="\t")
  ## The following test uses more than 2GB of memory, which is more
  ## than what 32-bit Windows can handle:
  is_32bit_windows <- .Platform$OS.type == "windows" &&
    .Platform$r_arch == "i386"
  if (!is_32bit_windows) {
    #Run MungeSumstats code
    reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter=FALSE,
                                                  bi_allelic_filter=FALSE,
                                                  allele_flip_check=FALSE)
    res_dt <- data.table::fread(reformatted)
    #correct names of MungeSumstats::eduAttainOkbay
    names(sumstats_dt) <- c("SNP","CHR","BP","A1","A2","FRQ","Beta","SE","P")
    #get order same
    setkey(res_dt,SNP)
    setkey(sumstats_dt,SNP)
    #add CHR, BP to org
    sumstats_dt[res_dt,CHR_der:=i.CHR]
    sumstats_dt[res_dt,BP_der:=i.BP]
    #remove any that weren't found in reference
    sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt),]
    expect_equal(sumstats_dt$BP,sumstats_dt$BP_der)
    expect_equal(sumstats_dt$CHR,sumstats_dt$CHR_der)
  }
  else{
    expect_equal(is_32bit_windows,TRUE)
    expect_equal(is_32bit_windows,TRUE)
  }
})
