test_that("SNPs not on reference genome are removed", {
    file <- tempfile()
    # Update ID from line 3 to check it is deleted -
    # "rs79925071" is not on ref genome GRCh37
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    eduAttainOkbay_missing <- eduAttainOkbay
    eduAttainOkbay_missing[3] <-
        "rs79925071\t6\t98584733\tA\tC\t0.5019\t0.024\t0.003\t2.457e-19"
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay_missing, con = file)
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" &&
        .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = TRUE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            log_folder_ind = TRUE
        )
        reformatted_lines <- readLines(reformatted$sumstats)
        # Should equal org apart from this one line
        writeLines(eduAttainOkbay, con = file)
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = TRUE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE
        )
        org_lines <- readLines(org)
        # reordering in function, line 3 rs9320913 is now 58
        # expect_equal(setequal(reformatted_lines,org_lines[-58]),TRUE)
        expect_equal(setequal(reformatted_lines, org_lines), TRUE)
        
        #also check get genome builds works
        eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
             package = "MungeSumstats")
        sumstats_list <- list(ss1 = eduAttainOkbayPth)
        ref_genomes <- get_genome_builds(sumstats_list = sumstats_list,
                                         sampled_snps = 50)
        expect_equal(all.equal(ref_genomes,list("ss1"="GRCH37")),TRUE)
        
    } else {
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})
