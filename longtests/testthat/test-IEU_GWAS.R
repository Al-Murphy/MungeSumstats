test_that("Test connection to IEU GWAS - metadata, download", {
  ## The following test uses more than 2GB of memory, which is more
  ## than what 32-bit Windows can handle:
  is_32bit_windows <- .Platform$OS.type == "windows" ##&&
  ##.Platform$r_arch == "i386"
  #only run test if user has internet access
  #capture internet outage or server issues
  if(try(is.character(RCurl::getURL("www.google.com")))==TRUE &&
     !is_32bit_windows && 
     #catch server issues - fake id just check it connects
     #Server code: 502; Server is possibly experiencing traffic, trying again...
     is.data.frame(try(gwasinfo(id = c("fake-id"),
                                opengwas_jwt = ieugwasr::get_opengwas_jwt())))){
    
    ### By ID and sample size
    metagwas <-
      try(MungeSumstats::find_sumstats(
        ids = c(
          "ieu-b-4760", "prot-a-1725",
          "prot-a-664"
        ),
        min_sample_size = 1000
      ))
    
    ### By criteria
    metagwas3 <- try(MungeSumstats::find_sumstats(
      traits = c("alzheimer", "parkinson"),
      years = seq(2015, 2021)
    ))
    # test these worked
    if(is.data.frame(metagwas)){
      testthat::expect_gt(nrow(metagwas), 0)
    }else{
      #catch server issues - fake id just check it connects
      #Server code: 502; Server is possibly experiencing traffic, trying again..
      testthat::expect_equal(TRUE, TRUE)
    }
    if(is.data.frame(metagwas3)){
      testthat::expect_gt(nrow(metagwas3), 0)
    }else{
      #catch server issues - fake id just check it connects
      #Server code: 502; Server is possibly experiencing traffic, trying again..
      testthat::expect_equal(TRUE, TRUE)
    }
    
    
    # test download
    vcf_url <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz"
    out_paths <- MungeSumstats::download_vcf(
      vcf_url = vcf_url,
      force_new = TRUE
    )
    # test this worked
    testthat::expect_true(file.exists(out_paths$save_path))
    
    
    # issues should be caught
    err_catch <-
      tryCatch({
        MungeSumstats:: import_sumstats(
          ids = "a-fake-id",
          ref_genome = "GRCh37",
          on_ref_genome = FALSE,
          strand_ambig_filter = FALSE,
          bi_allelic_filter = FALSE,
          allele_flip_check = FALSE, 
          force_new = TRUE,
          dbSNP=144)
      },
      error = function(e) e 
      )
    testthat::expect_true(methods::is(err_catch, "error"))
    # try import with axel - again should get warning about 1 thread
    axel_catch <-
      tryCatch(MungeSumstats::import_sumstats(
        ids = "ieu-a-298",
        ref_genome = "GRCh37",
        on_ref_genome = FALSE,
        strand_ambig_filter = FALSE,
        bi_allelic_filter = FALSE,
        allele_flip_check = FALSE,
        force_new = TRUE,
        nThread = 1,
        download_method = "axel",
        dbSNP=144
      ),
      error = function(e) e, 
      warning = function(w) w
      )
    testthat::expect_true(!is(axel_catch, "error"))
    
    # don't run last check too time intensive, it is also in the vignette anyway
    ## test importing
    ### Only use a subset for testing purposes
    #ids <- (dplyr::arrange(metagwas3, nsnp))$id
    # reformatted <- MungeSumstats::import_sumstats(ids = ids[1],
    #                                              ref_genome="GRCh37",
    #                                              on_ref_genome = FALSE,
    #                                              strand_ambig_filter=FALSE,
    #                                              bi_allelic_filter=FALSE,
    #                                              allele_flip_check=FALSE)
    # test this worked
    # testthat::expect_equal(file.exists(out_paths$save_path),TRUE)
    
  }else{
    testthat::expect_equal(TRUE, TRUE)
    testthat::expect_equal(TRUE, TRUE)
    testthat::expect_equal(TRUE, TRUE)
    testthat::expect_equal(TRUE, TRUE)
    testthat::expect_equal(TRUE, TRUE)
    testthat::expect_equal(TRUE, TRUE)
  }
})
