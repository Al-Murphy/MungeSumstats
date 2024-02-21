test_that("Handle dup cols", {
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
        # read it in and make N
        sumstats_dt <- data.table::fread(file)
        # Add dup column and make it not an integer
        sumstats_dt[, P := Pval]
        data.table::fwrite(x = sumstats_dt, file = file, sep = "\t")
        # rename Pval to p
        sumstats_txt <- readLines(file)
        sumstats_txt[1] <-
            "MarkerName\tCHR\tPOS\tA1\tA2\tEAF\tBeta\tSE\tP\tP"
        writeLines(sumstats_txt, file)
        # Run MungeSumstats code
        # Don't convert n to integers as this may round down
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            N_std = 5, convert_n_int = FALSE,
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            log_folder_ind = TRUE,
            imputation_ind = TRUE,
            log_mungesumstats_msgs = TRUE,
            dbSNP=144
        )
    
        # first check log file
        expect_equal(file.exists(reformatted$log_files$MungeSumstats_log_msg), 
                        TRUE)
        expect_equal(
            file.exists(reformatted$log_files$MungeSumstats_log_output),
            TRUE
        )
    
        res_dt <- data.table::fread(reformatted$sumstats)
    
        # should equal org
        writeLines(eduAttainOkbay, con = file)
        org_before_dup <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            N_std = 5, convert_n_int = FALSE,
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            log_folder_ind = TRUE,
            imputation_ind = TRUE,
            dbSNP=144
        )
        org_dt <- data.table::fread(org_before_dup$sumstats)
        expect_equal(all.equal(org_dt, res_dt), TRUE)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})
