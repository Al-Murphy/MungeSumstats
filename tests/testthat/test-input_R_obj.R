test_that("Can pass in R objects of summary statistics", {
    file <- tempfile()
    # write the Educational Attainment GWAS to a temp file for testing
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    writeLines(eduAttainOkbay, con = file)
    # read it in and combine CHR BP columns
    sumstats_dt <- data.table::fread(file)
    # Run MungeSumstats code
    path_return <- MungeSumstats::format_sumstats(file,
        ref_genome = "GRCh37",
        on_ref_genome = FALSE,
        strand_ambig_filter = FALSE,
        bi_allelic_filter = FALSE,
        allele_flip_check = FALSE
    )
    dt_return <- MungeSumstats::format_sumstats(sumstats_dt,
        ref_genome = "GRCh37",
        on_ref_genome = FALSE,
        strand_ambig_filter = FALSE,
        bi_allelic_filter = FALSE,
        allele_flip_check = FALSE
    )
    sumstats_rtrn_path <- data.table::fread(path_return)
    sumstats_rtrn_dt <- data.table::fread(dt_return)
    expect_equal(all.equal(sumstats_rtrn_dt, sumstats_rtrn_path), TRUE)
})
