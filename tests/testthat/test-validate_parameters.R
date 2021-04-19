test_that("validate parameters", {
  
  file <- tempfile()
  #write the Educational Attainment GWAS to a temp file for testing
  writeLines(MungeSumstats::eduAttainOkbay,con = file)
  
  #Check that submitting incorrect parameters gives errors
  #Run MungeSumstats code
  #passing incorrect path
  not_file <- tempfile()
  error_return <-
    tryCatch( MungeSumstats::format_sumstats(not_file,ref_genome="GRCh37",
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  #incorrect ref genome
  error_return2 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="fake",
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  #incorrect convert_small_p value
  error_return3 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                              convert_small_p = "fake",
                                              on_ref_genome = FALSE,
                                              strand_ambig_filter=FALSE,
                                              bi_allelic_filter=FALSE,
                                              allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  #incorrect convert_n_int value
  error_return4 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             convert_n_int = "fake",
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect on_ref_genome value
  error_return5 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             on_ref_genome = "YES",
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect strand_ambig_filter value
  error_return6 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter="YES",
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect INFO_filter value
  error_return7 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             INFO_filter=TRUE,
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect N_std value
  error_return8 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             N_std=FALSE,
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect rmv_chr value
  error_return9 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             rmv_chr=c("X","FAKE"),
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter=FALSE,
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect bi_allelic value
  error_return10 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             bi_allelic_filter="FAKE",
                                             allele_flip_check=FALSE),
              error = function(e) e,
              warning = function(w) w
    )
  
  #incorrect allele_flip_check value
  error_return11 <-
    tryCatch( MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                             on_ref_genome = FALSE,
                                             strand_ambig_filter=FALSE,
                                             allele_flip_check="FAKE"),
              error = function(e) e,
              warning = function(w) w
    )
  
  expect_equal(is(error_return, "error"), TRUE)
  expect_equal(is(error_return2, "error"), TRUE)
  expect_equal(is(error_return3, "error"), TRUE)
  expect_equal(is(error_return4, "error"), TRUE)
  expect_equal(is(error_return5, "error"), TRUE)
  expect_equal(is(error_return6, "error"), TRUE)
  expect_equal(is(error_return7, "error"), TRUE)
  expect_equal(is(error_return8, "error"), TRUE)
  expect_equal(is(error_return9, "error"), TRUE)
  expect_equal(is(error_return10, "error"), TRUE)
  expect_equal(is(error_return11, "error"), TRUE)
})
