test_that("Test that user inputted mapping file is appropriate", {
    file <- tempfile()
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay, con = file)

    # define mapping file - leave out some essential columns - A2
    data("sumstatsColHeaders")
    essential_cols <- c("SNP", "CHR", "BP", "P", "A1", "A2")
    signed_cols <- c("Z", "OR", "BETA", "LOG_ODDS", "SIGNED_SUMSTAT")
    user_map <- sumstatsColHeaders[sumstatsColHeaders$Corrected %in%
        essential_cols[-6] |
        sumstatsColHeaders$Corrected %in%
            signed_cols, ]
    # Run MungeSumstats code
    fail_return <-
        tryCatch(MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            mapping_file = user_map
        ),
        error = function(e) e,
        warning = function(w) w
        )
    expect_equal(is(fail_return, "error"), TRUE)
    # update user map to contain all essential columns
    user_map2 <- sumstatsColHeaders[sumstatsColHeaders$Corrected %in%
        essential_cols |
        sumstatsColHeaders$Corrected %in%
            signed_cols, ]
    # This time should run
    reformatted <- MungeSumstats::format_sumstats(file,
        ref_genome = "GRCh37",
        on_ref_genome = FALSE,
        strand_ambig_filter = FALSE,
        bi_allelic_filter = FALSE,
        allele_flip_check = FALSE,
        mapping_file = user_map2,
        return_data = TRUE,
        return_format = "GRanges"
    )


    expect_equal(is(reformatted, "GRanges"), TRUE)
})
