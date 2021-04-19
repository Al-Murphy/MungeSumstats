test_that("Handle small p-values", {
  file <- tempfile()
  #Remove data from line 3 to check it is deleted
  eduAttainOkbay_missing <- MungeSumstats::eduAttainOkbay
  eduAttainOkbay_missing[3] <-
    "rs9320913\t6\t98584733\tA\tC\t0.5019\t0.024\t0.003\t2.457e-339" 
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(eduAttainOkbay_missing,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                on_ref_genome = FALSE,
                                                strand_ambig_filter=FALSE,
                                                bi_allelic_filter=FALSE,
                                                allele_flip_check=FALSE)
  reformatted_lines <- readLines(reformatted)
  #Should equal org apart from this one line
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  org <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                        on_ref_genome = FALSE,
                                        strand_ambig_filter=FALSE,
                                        bi_allelic_filter=FALSE,
                                        allele_flip_check=FALSE)
  org_lines <- readLines(org)
  #rows get reordered in function so 3 -> 58
  expect_equal(reformatted_lines[-58],org_lines[-58])
})
