test_that("Handle n when its 5 std dev > mean", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  eduAttainOkbay <- readLines(system.file("extdata","eduAttainOkbay.txt",
                                          package="MungeSumstats"))
  writeLines(eduAttainOkbay,con = file)
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
                                                  allele_flip_check=FALSE,
                                                  log_folder_ind = TRUE)
  res_dt <- data.table::fread(reformatted$sumstats)
  expect_equal(!(rmv_snp %in% res_dt$SNP),TRUE)
  
  #Run MungeSumstats code
  #Don't convert n to integers as this may round down
  #set N_dropNA to TRUE
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                N_std=5,convert_n_int=FALSE,
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE,
                                                log_folder_ind = TRUE,
                                                N_dropNA = TRUE)
  res_dt <- data.table::fread(reformatted$sumstats)
  
  #set N_dropNA to FALSE
  reformatted2 <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                N_std=5,convert_n_int=FALSE,
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE,
                                                log_folder_ind = TRUE,
                                                N_dropNA = FALSE)
  res_dt2 <- data.table::fread(reformatted2$sumstats)
  expect_equal(!(rmv_snp %in% res_dt$SNP),TRUE)
  #should be no na's
  expect_equal(all.equal(res_dt,res_dt2),TRUE)
})
