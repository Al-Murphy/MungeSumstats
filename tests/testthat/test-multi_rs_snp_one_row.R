test_that("Handle more than 1 rs IDs in one row", {
  ## The following test uses more than 2GB of memory, which is more
  ## than what 32-bit Windows can handle:
  is_32bit_windows <- .Platform$OS.type == "windows" #&&
  #.Platform$r_arch == "i386"
  if (!is_32bit_windows  && Sys.info()["sysname"]=="Linux") {
    file <- tempfile()
    # Remove data from line 3 to check it is deleted
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
                                            package = "MungeSumstats"
    ))
    eduAttainOkbay_missing <- eduAttainOkbay
    eduAttainOkbay_missing[3] <-
      "rs9320913_rs1234_rs_45678\t6\t98584733\tA\tC\t0.5019\t0.024\t0.003\t2.457e-19"
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay_missing, con = file)
    
    
    # make changes for log check
    file2 <- tempfile()
    # multiple rs id already there
    # data already has 8 strand ambiguous snps, 1 snp A1 A2 doesn't match ref gen
    # add missing data(EAF)
    eduAttainOkbay_missing[2] <-
      "rs12987662\t2\t100821548\tA\tC\t\t0.027\t0.003\t2.693e-24"
    # make beta 0
    eduAttainOkbay_missing[4] <-
      "rs11712056\t3\t49914397\tT\tC\t0.5504\t0\t0.003\t3.304e-19"
    # make se negative
    eduAttainOkbay_missing[9] <-
      "rs2456973\t12\t56416928\tA\tC\t0.6791\t-0.02\t-0.003\t1.064e-12"
    # make duplicate SNP IDs
    eduAttainOkbay_missing[10] <-
      "rs165633\t12\t123767929\tA\tG\t0.2257\t0.023\t0.003\t1.258e-12"
    # make duplicate base pair positions - bp from row 14
    eduAttainOkbay_missing[13] <-
      "rs11191193\t4\t140764124\tA\tG\t0.6511\t0.018\t0.003\t5.444e-11"
    # make snp missing rs but change chr and bp to not real ones so then removed
    eduAttainOkbay_missing[15] <-
      "11210860\t1\t-1\tA\tG\t0.3694\t0.017\t0.003\t2.359e-10"
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay_missing, con = file2)
    # Run MungeSumstats code
    reformatted <-
      MungeSumstats::format_sumstats(file,
                                     ref_genome = "GRCh37",
                                     on_ref_genome = FALSE,
                                     strand_ambig_filter = FALSE,
                                     bi_allelic_filter = FALSE,
                                     allele_flip_check = FALSE,
                                     imputation_ind = TRUE,
                                     remove_multi_rs_snp = FALSE,
                                     dbSNP=144
      )
    reformatted_lines <- data.table::fread(reformatted)
    # Should equal org apart from this one line
    writeLines(eduAttainOkbay, con = file)
    org <- MungeSumstats::format_sumstats(file,
                                          ref_genome = "GRCh37",
                                          on_ref_genome = FALSE,
                                          strand_ambig_filter = FALSE,
                                          bi_allelic_filter = FALSE,
                                          allele_flip_check = FALSE,
                                          imputation_ind = TRUE,
                                          remove_multi_rs_snp = FALSE,
                                          dbSNP=144
    )
    org_lines <- data.table::fread(org)
    
    # remove imputation column
    reformatted_lines[, convert_multi_rs_SNP := NULL]
    # reordering makes line 3 got to 58
    testthat::expect_equal(reformatted_lines, org_lines)
    
    # check log files
    # Run MungeSumstats code
    reformatted_log <-
      MungeSumstats::format_sumstats(file2,
                                     ref_genome = "GRCh37",
                                     on_ref_genome = FALSE,#TRUE,
                                     strand_ambig_filter = TRUE,
                                     bi_allelic_filter = FALSE,#TRUE,
                                     allele_flip_check = FALSE,#TRUE,
                                     imputation_ind = TRUE,
                                     remove_multi_rs_snp = TRUE,
                                     log_folder_ind = TRUE,
                                     dbSNP=144
      )
    
    # expect 5 log files
    testthat::expect_equal(length(reformatted_log$log_files), 5)
    # next check number of rows in each
    results <- c()
    for (log_i in reformatted_log$log_files) {
      data_log_i <- data.table::fread(log_i)
      if (grepl("snp_strand_ambiguous", log_i)) {
        results <- c(results, nrow(data_log_i) == 8)
      } else {
        results <- c(results, nrow(data_log_i) == 1|
                       nrow(data_log_i) == 0)
      }
    }
    expect_equal(all(results), TRUE)
  } else {
    expect_equal((is_32bit_windows||!Sys.info()["sysname"]=="Linux"), TRUE)
    expect_equal((is_32bit_windows||!Sys.info()["sysname"]=="Linux"), TRUE)
    expect_equal((is_32bit_windows||!Sys.info()["sysname"]=="Linux"), TRUE)
  }
})