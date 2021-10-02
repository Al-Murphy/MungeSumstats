test_that("will fail without a vital column", {
    file <- tempfile()
    # write the Educational Attainment GWAS to a temp file for testing
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    writeLines(eduAttainOkbay, con = file)
    # read it in and combine CHR BP columns
    sumstats_dt <- data.table::fread(file, nThread = 1)
    # Remove vital column - P
    sumstats_dt[, Pval := NULL]
    data.table::fwrite(
        x = sumstats_dt,
        file = file,
        sep = "\t",
        nThread = 1
    )
    # Run MungeSumstats code
    error_return <-
        tryCatch(MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE
        ),
        error = function(e) e,
        warning = function(w) w
        )
    expect_equal(is(error_return, "error"), TRUE)
})
