test_that("Test connection to IEU GWAS - metadata, download", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
    ##.Platform$r_arch == "i386"
    #only run test if user has internet access
    #capture internet outage or server issues
    if(try(is.character(getURL("www.google.com")))==TRUE &&
       !is_32bit_windows){
        ### By ID
        metagwas <-
            MungeSumstats::find_sumstats(ids = c(
                "ieu-b-4760", "prot-a-1725",
                "prot-a-664"
            ))
    
        ### By ID amd sample size
        metagwas2 <-
            MungeSumstats::find_sumstats(
                ids = c(
                    "ieu-b-4760", "prot-a-1725",
                    "prot-a-664"
                ),
                min_sample_size = 1000
            )
    
        ### By criteria
        metagwas3 <- MungeSumstats::find_sumstats(
            traits = c("alzheimer", "parkinson"),
            years = seq(2015, 2021)
        )
        # test these worked
        testthat::expect_gt(nrow(metagwas), 0)
        testthat::expect_gt(nrow(metagwas2), 0)
        testthat::expect_gt(nrow(metagwas3), 0)
    
    
        # test download
        vcf_url <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz"
        out_paths <- MungeSumstats::download_vcf(
            vcf_url = vcf_url,
            force_new = TRUE
        )
        # test this worked
        testthat::expect_true(file.exists(out_paths$save_path))
    
        # test importing
        ### Only use a subset for testing purposes
        ids <- (dplyr::arrange(metagwas3, nsnp))$id
        # issues should be caught
        err_catch <-
            tryCatch(MungeSumstats::import_sumstats(
                ids = "a-fake-id",
                ref_genome = "GRCh37",
                on_ref_genome = FALSE,
                strand_ambig_filter = FALSE,
                bi_allelic_filter = FALSE,
                allele_flip_check = FALSE,
                force_new = TRUE
            ),
            error = function(e) e, warning = function(w) w
            )
        testthat::expect_true(is(err_catch, "error"))
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
                download_method = "axel"
            ),
            error = function(e) e, warning = function(w) w
            )
        testthat::expect_true(!is(axel_catch, "error"))

        # don't run last check too time intensive, it is also in the vignette anyway
        # reformatted <- MungeSumstats::import_sumstats(ids = ids[1],
        #                                              ref_genome="GRCh37",
        #                                              on_ref_genome = FALSE,
        #                                              strand_ambig_filter=FALSE,
        #                                              bi_allelic_filter=FALSE,
        #                                              allele_flip_check=FALSE)
        # test this worked
        # testthat::expect_equal(file.exists(out_paths$save_path),TRUE)
    }
    else{
        expect_equal(TRUE, TRUE)
    }
})
