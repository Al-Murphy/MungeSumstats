test_that("Can correctly separate CHR:BP:A2:A1, two columns", {
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
        # read it in and combine CHR BP columns
        sumstats_dt <- data.table::fread(file)
        # Keep Org to validate values
        sumstats_dt_missing <- data.table::copy(sumstats_dt)
        sumstats_dt_missing[, CHR_BP_A2_A1 := paste0(CHR, ":", POS, ":", 
                                                        A2, ":", A1)]
        sumstats_dt_missing[, CHR := NULL]
        sumstats_dt_missing[, POS := NULL]
        sumstats_dt_missing[, A1 := NULL]
        sumstats_dt_missing[, A2 := NULL]
        data.table::fwrite(x = sumstats_dt_missing, file = file, sep = "\t")
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        res_dt <- data.table::fread(reformatted)
        # Should give same result as separated
        data.table::fwrite(x = sumstats_dt, file = file, sep = "\t")
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        org_dt <- data.table::fread(org)
        expect_equal(all.equal(org_dt, res_dt, ignore.row.order = TRUE,
                               ignore.col.order = TRUE), TRUE)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})
