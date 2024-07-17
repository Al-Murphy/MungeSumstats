test_that("log10 p-value", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
      # read it in and make N
      sumstats_dt <- data.table::fread(system.file("extdata", "eduAttainOkbay.txt",
                                                   package = "MungeSumstats"), 
                                       nThread = 1)
      #test -log10
      sumstats_dt[,LP:=-1*log(Pval,base=10)]
      data.table::setnames(sumstats_dt,"Pval","Pval_org")
      reformatted <- MungeSumstats::format_sumstats(sumstats_dt,
                                                    ref_genome = "GRCh37",
                                                    INFO_filter = 0.9,
                                                    on_ref_genome = FALSE,
                                                    strand_ambig_filter = FALSE,
                                                    bi_allelic_filter = FALSE,
                                                    allele_flip_check = FALSE,
                                                    log_folder_ind = FALSE,
                                                    imputation_ind = FALSE,
                                                    dbSNP=144,
                                                    ignore_multi_trait=TRUE,
                                                    return_data = TRUE)
      testthat::expect_equal(reformatted$PVAL_ORG,reformatted$P)
      #test log10
      sumstats_dt <- data.table::fread(system.file("extdata", "eduAttainOkbay.txt",
                                                   package = "MungeSumstats"), 
                                       nThread = 1)
      sumstats_dt[,LP:=log(Pval,base=10)]
      data.table::setnames(sumstats_dt,"Pval","Pval_org")
      reformatted <- MungeSumstats::format_sumstats(sumstats_dt,
                                                    ref_genome = "GRCh37",
                                                    INFO_filter = 0.9,
                                                    on_ref_genome = FALSE,
                                                    strand_ambig_filter = FALSE,
                                                    bi_allelic_filter = FALSE,
                                                    allele_flip_check = FALSE,
                                                    log_folder_ind = FALSE,
                                                    imputation_ind = FALSE,
                                                    dbSNP=144,
                                                    ignore_multi_trait=TRUE,
                                                    return_data = TRUE)
      testthat::expect_equal(reformatted$PVAL_ORG,reformatted$P)
    }    
    else{
      testthat::expect_equal(is_32bit_windows, TRUE)
      testthat::expect_equal(is_32bit_windows, TRUE)
    }
})
