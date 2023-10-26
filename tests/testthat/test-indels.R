test_that("non-biallelic SNPs are removed", {
  ## The following test uses more than 2GB of memory, which is more
  ## than what 32-bit Windows can handle:
  is_32bit_windows <- .Platform$OS.type == "windows" #&&
  #.Platform$r_arch == "i386"
  if (!is_32bit_windows && Sys.info()["sysname"]=="Linux") {
    #test to ensure indels aren't removed
    # also test indel missing RS ID removed rather than imputing wrong RS ID
    ss_indel <- data.table::data.table("SNP"=c("rs34589910","rs12987662",
                                               "4:6364621"),
                                       "CHR"=c(4,2,4),
                                       "BP"=c(6364621,100821548,6364621),
                                       "A1"=c("C","A","C"),
                                       "A2"=c("CG","C","CG"),
                                       "Uniq.a1a2"=c("4:6364621_C_CG","aa",
                                                     "4:6364621_C_CG"),
                                       "EAF"=c(0.0945334,0.3787,0.0945334),
                                       "BETA"=c(-0.00625732297153778,0.027,
                                                -0.00625732297153778),
                                       "P"=c(0.4883341,2.693e-24,0.4883341))
    
    reformatted_ss_ad <-
      MungeSumstats::format_sumstats(ss_indel,ref_genome="GRCh37",
                                     convert_small_p=TRUE,
                                     allele_flip_check=TRUE,
                                     snp_ids_are_rs_ids=TRUE,
                                     return_data=TRUE,
                                     nThread=2,
                                     on_ref_genome = TRUE,
                                     indels = TRUE,
                                     log_folder_ind = TRUE,
                                     dbSNP=144)
    #SNP ID is an indel so won't exist in our SNP reference dataset
    testthat::expect_equal("rs34589910" %in% 
                             reformatted_ss_ad$sumstats$SNP,TRUE)
    #check that indel missing RS ID is removed rather than imputing 
    testthat::expect_equal(nrow(fread(
      reformatted_ss_ad$log_files$snp_missing_rs)),1)
    
  } else {
    testthat::expect_equal((is_32bit_windows||
                              !Sys.info()["sysname"]=="Linux"), TRUE)
    testthat::expect_equal((is_32bit_windows||
                              !Sys.info()["sysname"]=="Linux"), TRUE)
  }
})
