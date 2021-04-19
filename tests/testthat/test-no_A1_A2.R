test_that("Imputation of A1/A2 correctly", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and drop CHR BP columns
  sumstats_dt <- data.table::fread(file)
  #Keep Org to validate values
  sumstats_dt_missing <- data.table::copy(sumstats_dt)
  sumstats_dt_missing[,A1:=NULL]
  sumstats_dt_missing[,A2:=NULL]
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
    
    #Check with just one of A1 A2
    #write the Educational Attainment GWAS to a temp file for testing
    writeLines(MungeSumstats::eduAttainOkbay,con = file)
    #read it in and drop CHR BP columns
    sumstats_dt <- data.table::fread(file)
    #Keep Org to validate values
    sumstats_dt_missing <- data.table::copy(sumstats_dt)
    sumstats_dt_missing[,A1:=NULL]
    data.table::fwrite(x=sumstats_dt_missing, file=file, sep="\t")
    #Run MungeSumstats code
    reformatted2 <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter=FALSE,
                                                  bi_allelic_filter=FALSE,
                                                  allele_flip_check=FALSE)
    res_dt2 <- data.table::fread(reformatted2)
    
    #correct names of MungeSumstats::eduAttainOkbay
    names(sumstats_dt) <- c("SNP","CHR","BP","A1","A2","FRQ","Beta","SE","P")
    #get order same
    setkey(res_dt,SNP)
    setkey(res_dt2,SNP)
    setkey(sumstats_dt,SNP)
    #add A1 to org
    sumstats_dt[res_dt,A1_der:=i.A1]
    #add A2 to org
    sumstats_dt[res_dt,A2_der:=i.A2]
    
    #add second A2 to org
    sumstats_dt[res_dt2,A2_der2:=i.A2]
    
    #remove any that weren't found in reference
    sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt),]
    #random chance would be 25% - so this is significantly greater than that
    A1_valid <- mean(sumstats_dt$A1==sumstats_dt$A1_der)>0.5 #50% threshold
    #expect A1s to be equal, won't be always right
    expect_equal(A1_valid,TRUE)
    A2_valid <- mean(sumstats_dt$A2==sumstats_dt$A2_der)>0.45 #45% threshold
    #expect A2s to be equal, all won't be since can be multiple A2s
    expect_equal(A2_valid,TRUE)
    
    A2_valid2 <- mean(sumstats_dt$A2==sumstats_dt$A2_der2)>0.45 #45% threshold
    #expect A2s to be equal, all won't be since can be multiple A2s
    expect_equal(A2_valid2,TRUE)
  }
  else{
    expect_equal(is_32bit_windows,TRUE)
    expect_equal(is_32bit_windows,TRUE)
    expect_equal(is_32bit_windows,TRUE)
  }
})
