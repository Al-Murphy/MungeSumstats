test_that("VCF is correctly formatted", {
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::ieuAmlVcf,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE)
  reformatted_lines <- readLines(reformatted)
  #check manually against first five SNPs
  corr_res <- c(
    "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tLP\tAF\tP",
    "rs58108140\t1\t10583\tG\tA\t0.0312\t0.0393\t0.369267\t0.1589\t0.427300105456596",
    "rs806731\t1\t30923\tG\tT\t-0.0114\t0.0353\t0.126854\t0.7843\t0.746699739815279",
    "rs116400033\t1\t51479\tT\tA\t0.0711\t0.037\t1.26241\t0.1829\t0.0546499790752282",
    "rs146477069\t1\t54421\tA\tG\t-0.024\t0.083\t0.112102\t0.0352\t0.772499131799648")

  expect_equal(reformatted_lines[1:5],corr_res)
})
