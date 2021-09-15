test_that("Filter SNPs where INFO<0.9", {
    file <- tempfile()
    # write the Educational Attainment GWAS to a temp file for testing
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    writeLines(eduAttainOkbay, con = file)
    # read it in and make N
    sumstats_dt <- data.table::fread(file, nThread = 1)
    # Add N column and make it not an integer
    set.seed(101)
    sumstats_dt[, INFO := runif(nrow(sumstats_dt)) * 2]
    # get SNPs with INFO<0.9
    rmv_snps <- sumstats_dt[INFO < 0.9, ]$MarkerName
    data.table::fwrite(x = sumstats_dt, file = file, sep = "\t")
    # Run MungeSumstats code
    reformatted <- MungeSumstats::format_sumstats(file,
        ref_genome = "GRCh37",
        INFO_filter = 0.9,
        on_ref_genome = FALSE,
        strand_ambig_filter = FALSE,
        bi_allelic_filter = FALSE,
        allele_flip_check = FALSE,
        log_folder_ind = TRUE,
        imputation_ind = TRUE
    )
    res_dt <- data.table::fread(reformatted$sumstats, nThread = 1)
    testthat::expect_equal(all(!rmv_snps %in% res_dt$SNP), TRUE)
})
