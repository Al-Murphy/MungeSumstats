test_that("Test that if a0/a1 passed, formatted correctly", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # The dataset's alleles need to be flipped to test
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
            package = "MungeSumstats"
        ))
        eduAttainOkbay_missing <- eduAttainOkbay
        # A0/A1 go to ref/alt whereas A1/A2 goes to ref/alt
        eduAttainOkbay_missing[1] <-
            "MarkerName\tCHR\tPOS\tA0\tA1\tEAF\tBeta\tSE\tPval"
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay_missing, con = file)
    
    
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE
        )
        reformatted_lines <- readLines(reformatted)
        # Should equal org since the effect should be corrected
        writeLines(eduAttainOkbay, con = file)
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE
        )
        org_lines <- readLines(org)
        # reordering in function
        expect_equal(setequal(reformatted_lines, org_lines), TRUE)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})
