test_that("Can handle space delimited files", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # write the Educational Attainment GWAS to a temp file for testing
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
            package = "MungeSumstats"
        ))
        writeLines(eduAttainOkbay, con = file)
        # read it in and write as space delimited
        sumstats_dt <- data.table::fread(file, nThread = 1)
        data.table::fwrite(x = sumstats_dt, file = file, sep = " ")
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        res_dt <- data.table::fread(reformatted, nThread = 1)
        # check against results of normal run should be the exact same
        file2 <- tempfile()
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay, con = file2)
        org <- MungeSumstats::format_sumstats(file2,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        org_dt <- data.table::fread(org, nThread = 1)
        expect_equal(res_dt, org_dt)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})
