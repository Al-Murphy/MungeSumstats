test_that("Test that allele columns and effect columns flipped correctly", {
    file <- tempfile()
    # The dataset's alleles need to be flipped to test
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    eduAttainOkbay_missing <- eduAttainOkbay
    eduAttainOkbay_missing[1] <-
        "MarkerName\tCHR\tPOS\tA2\tA1\tEAF\tBeta\tSE\tPval"
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay_missing, con = file)
    # read in and manually change effect columns
    eduAttainOkbay_missing_dt <- data.table::fread(file)
    eduAttainOkbay_missing_dt[, Beta := Beta * -1]
    eduAttainOkbay_missing_dt[, EAF := 1 - EAF]
    data.table::fwrite(eduAttainOkbay_missing_dt,
        file = file, sep = "\t"
    )
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
        ##.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = TRUE,
            allele_flip_check = TRUE
        )
        reformatted_lines <- readLines(reformatted)
        # Should equal org since the effect should be corrected
        writeLines(eduAttainOkbay, con = file)
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = TRUE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = TRUE,
            allele_flip_check = TRUE
        )
        org_lines <- readLines(org)
        # reordering in function
        expect_equal(setequal(reformatted_lines, org_lines), TRUE)
    } else {
        expect_equal(is_32bit_windows, TRUE)
    }
})
