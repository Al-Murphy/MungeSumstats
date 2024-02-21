test_that("Handle small p-values", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # Remove data from line 3 to check it is deleted
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
            package = "MungeSumstats"
        ))
        eduAttainOkbay_missing <- eduAttainOkbay
        eduAttainOkbay_missing[3] <-
            "rs9320913\t6\t98584733\tA\tC\t0.5019\t0.024\t0.003\t2.457e-339"
        problem_snp <- "rs9320913"
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay_missing, con = file)
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        reformatted_lines <- readLines(reformatted)
        # Should equal org apart from this one line
        writeLines(eduAttainOkbay, con = file)
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        org_lines <- readLines(org)
        rsid_index <- grep(problem_snp, org_lines, ignore.case = TRUE)
        # rows get reordered in function so 3 -> 58
        expect_equal(reformatted_lines[-rsid_index], org_lines[-rsid_index])
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})
