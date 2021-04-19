test_that("Handle n when its 5 std dev > mean", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  #read it in and make N
  sumstats_dt <- data.table::fread(file)
  #Add N column and make it not an integer
  sumstats_dt[,N:=round(10*runif(nrow(sumstats_dt)),0)]
  #Ensure 1 value is > 5std dev above the mean
  #Remember creating and adding in this value will increase the mean and sd
  mean_N<- mean(sumstats_dt$N)
  sd_N<- stats::sd(sumstats_dt$N)
  bigger_value <- mean_N + (6*sd_N)
  data.table::set(sumstats_dt, i=1L, j="N", value=bigger_value)
  data.table::fwrite(x=sumstats_dt, file=file, sep="\t")
  #ensure this SNP is removed
  rmv_snp <- sumstats_dt$MarkerName[1]
  #Run MungeSumstats code
  #Don't convert n to integers as this may round down
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                  N_std=5,convert_n_int=FALSE,
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter=FALSE,
                                                  bi_allelic_filter=FALSE,
                                                  allele_flip_check=FALSE)
  res_dt <- data.table::fread(reformatted)
  expect_equal(!(rmv_snp %in% res_dt$SNP),TRUE)
})
