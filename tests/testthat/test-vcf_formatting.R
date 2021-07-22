test_that("VCF is correctly formatted", {
  file <- tempfile()
  #write the ALS GWAS, VCF file to a temp file for testing
  ALSvcf <- readLines(system.file("extdata","ALSvcf.vcf",
                                    package="MungeSumstats"))
  writeLines(ALSvcf,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE,
                                                INFO_filter = 0.01)
  reformatted_lines <- readLines(reformatted)
  #check manually against first five SNPs
  corr_res <- c(
    "SNP\tCHR\tBP\tA1\tA2\tINFO\tBETA\tSE\tLP\tAF\tP",
    "rs58108140\t1\t10583\tG\tA\t0.1589\t0.0312\t0.0393\t0.369267\t0.1589\t0.427300105456596",
    "rs806731\t1\t30923\tG\tT\t0.7843\t-0.0114\t0.0353\t0.126854\t0.7843\t0.746699739815279",
    "rs116400033\t1\t51479\tT\tA\t0.1829\t0.0711\t0.037\t1.26241\t0.1829\t0.0546499790752282",
    "rs146477069\t1\t54421\tA\tG\t0.0352\t-0.024\t0.083\t0.112102\t0.0352\t0.772499131799648")

  expect_equal(reformatted_lines[1:5],corr_res)
  
  #check allelic flipping with VCF
  file2 <- tempfile()
  #write the ALS GWAS, VCF file to a temp file for testing
  ALSvcf <- readLines(system.file("extdata","ALSvcf.vcf",
                                  package="MungeSumstats"))
  #update last SNP, flipping allelic direction
  snp_of_interest <- "rs146477069"
  rsid_index <- grep(snp_of_interest, ALSvcf, ignore.case = TRUE) 
  ALSvcf[rsid_index]<-
    "1\t54421\trs146477069\tG\tA\t.\tPASS\tAF=0.0352\tES:SE:LP:AF:ID\t+0.024:0.083:0.112102:0.0352:rs146477069"
  writeLines(ALSvcf,con = file2)
  
  #Run MungeSumstats code
  ## The following test uses more than 2GB of memory, which is more
  ## than what 32-bit Windows can handle:
  is_32bit_windows <- .Platform$OS.type == "windows" && 
    .Platform$r_arch == "i386"
  if (!is_32bit_windows) {
    reformatted_allelic_flip <- 
      MungeSumstats::format_sumstats(file2,ref_genome="GRCh37",
                                     on_ref_genome = FALSE,
                                     strand_ambig_filter=FALSE,
                                     bi_allelic_filter=FALSE,
                                     allele_flip_check=TRUE,
                                     allele_flip_drop = FALSE,
                                     INFO_filter = 0.01)
    reformatted_lines_af <- readLines(reformatted_allelic_flip)
    expect_equal(reformatted_lines,reformatted_lines_af)
    #also check outputing as different types
    pth <- system.file("extdata","ALSvcf.vcf",package="MungeSumstats")
    rtrn_dt <- MungeSumstats::format_sumstats(pth,ref_genome="GRCh37",
                                              on_ref_genome = FALSE,
                                              strand_ambig_filter=FALSE,
                                              bi_allelic_filter=FALSE,
                                              allele_flip_check=FALSE,
                                              allele_flip_drop = FALSE,
                                              INFO_filter = 0.01,
                                              return_data=TRUE,
                                              return_format="data.table")
    rtrn_grng <- MungeSumstats::format_sumstats(pth,ref_genome="GRCh37",
                                              on_ref_genome = FALSE,
                                              strand_ambig_filter=FALSE,
                                              bi_allelic_filter=FALSE,
                                              allele_flip_check=FALSE,
                                              allele_flip_drop = FALSE,
                                              INFO_filter = 0.01,
                                              return_data=TRUE,
                                              return_format="GRanges")
    rtrn_vrng <- MungeSumstats::format_sumstats(pth,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE,
                                                allele_flip_drop = FALSE,
                                                INFO_filter = 0.01,
                                                return_data=TRUE,
                                                return_format="VRanges")
    expect_equal(all(is(rtrn_grng)[1]=="GRanges",is(rtrn_vrng)[1]=="VRanges",
                     is(rtrn_dt)[1]=="data.table"),TRUE)
    
    #also test inferring the genome build
    rtrn_dt_infer <- MungeSumstats::format_sumstats(pth,
                                              on_ref_genome = FALSE,
                                              strand_ambig_filter=FALSE,
                                              bi_allelic_filter=FALSE,
                                              allele_flip_check=FALSE,
                                              allele_flip_drop = FALSE,
                                              INFO_filter = 0.01,
                                              return_data=TRUE,
                                              return_format="data.table")
    expect_equal(all.equal(rtrn_dt,rtrn_dt_infer),TRUE)
    
    #also test outputting ldsc_format ready format 
    rtrn_ldsc <- MungeSumstats::format_sumstats(pth,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE,
                                                allele_flip_drop = FALSE,
                                                INFO_filter = 0.01,
                                                ldsc_format = TRUE,
                                                compute_n = 1001)
    res<-data.table::fread(rtrn_ldsc)
    #check for necessary columns - https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
    ldsc_cols <- c("SNP","N","A1","A2","Z")
    expect_equal(all(ldsc_cols %in% names(res)),TRUE)
  }
  else{
    expect_equal(is_32bit_windows,TRUE)
    expect_equal(is_32bit_windows,TRUE)
    expect_equal(is_32bit_windows,TRUE)
    expect_equal(is_32bit_windows,TRUE)
  }  
  
  expect_equal(reformatted_lines[1:5],corr_res)
})
