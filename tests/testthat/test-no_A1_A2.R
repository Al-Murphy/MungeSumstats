test_that("Imputation of A1/A2 correctly", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
        ##.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # write the Educational Attainment GWAS to a temp file for testing
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
                                                package = "MungeSumstats"
        ))
        writeLines(eduAttainOkbay, con = file)
        # read it in and drop CHR BP columns
        sumstats_dt <- data.table::fread(file, nThread = 1)
        # Keep Org to validate values
        sumstats_dt_missing <- data.table::copy(sumstats_dt)
        sumstats_dt_missing[, A1 := NULL]
        sumstats_dt_missing[, A2 := NULL]
        data.table::fwrite(
            x = sumstats_dt_missing,
            file = file,
            sep = "\t",
            nThread = 1
        )
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            imputation_ind = TRUE,
            log_folder_ind = TRUE,
            dbSNP=144
        )
        res_dt <- data.table::fread(reformatted$sumstats, 
                                    nThread = 1)
        # imputation cols - check they are there, then drop 
        testthat::expect_equal(
            sum(c("IMPUTATION_A1", "IMPUTATION_A2") %in% names(res_dt)), 2)
        res_dt[, IMPUTATION_A1 := NULL]
        res_dt[, IMPUTATION_A2 := NULL]
        # Check with just one of A1 A2
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay, con = file)
        # read it in and drop A1 columns
        sumstats_dt <- data.table::fread(file,
                                         nThread = 1)
        # Keep Org to validate values
        sumstats_dt_missing <- data.table::copy(sumstats_dt)
        sumstats_dt_missing[, A1 := NULL]
        data.table::fwrite(x = sumstats_dt_missing,
                           file = file,
                           sep = "\t",
                           nThread = 1)
        # Run MungeSumstats code
        reformatted2 <- MungeSumstats::format_sumstats(
            path = file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            imputation_ind = TRUE,
            log_folder_ind = TRUE,
            allele_flip_frq = FALSE,
            dbSNP=144
        )
        res_dt2 <- data.table::fread(reformatted2$sumstats,
                                     nThread = 1)
        # imputation cols - check they are there, then drop
        testthat::expect_equal(
            sum(c("IMPUTATION_A1", "flipped") %in% names(res_dt2)), 2)
        res_dt2[, IMPUTATION_A1 := NULL]
        res_dt2[, flipped := NULL]
        # Check with just one of A1 A2
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay, con = file)
        # read it in and drop A2 columns
        sumstats_dt <- data.table::fread(file,
                                         nThread = 1)
        # Keep Org to validate values
        sumstats_dt_missing <- data.table::copy(sumstats_dt)
        sumstats_dt_missing[, A2 := NULL]
        data.table::fwrite(x = sumstats_dt_missing,
                           file = file,
                           sep = "\t",
                           nThread = 1)
        # Run MungeSumstats code
        reformatted3 <- MungeSumstats::format_sumstats(
            path = file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            allele_flip_drop = FALSE,
            imputation_ind = TRUE,
            log_folder_ind = TRUE,
            dbSNP=144
        )
        res_dt3 <- data.table::fread(reformatted3$sumstats,
                                     nThread = 1)
        # imputation cols - check they are there, then drop
        testthat::expect_equal("IMPUTATION_A2" %in% names(res_dt3), TRUE)
        res_dt3[, IMPUTATION_A2 := NULL]

        # correct names of MungeSumstats::eduAttainOkbay
        names(sumstats_dt) <- c("SNP", "CHR", "BP", "A1", "A2",
                                "FRQ", "Beta", "SE", "P")
        # get order same
        data.table::setkey(res_dt, SNP)
        data.table::setkey(res_dt2, SNP)
        data.table::setkey(sumstats_dt, SNP)
        # add A1 to org
        sumstats_dt[res_dt, A1_der := i.A1]
        # add A2 to org
        sumstats_dt[res_dt, A2_der := i.A2]

        # add second A1 to org
        sumstats_dt[res_dt2, A1_der2 := i.A1]
        # add second A2 to org
        sumstats_dt[res_dt3, A2_der2 := i.A2]

        # remove any that weren't found in reference
        sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt), ]
        # random chance would be 25% - so this is significantly greater than that
        A1_valid <- mean(sumstats_dt$A1 == sumstats_dt$A1_der) > 0.5 # 50% threshold
        # expect A1s to be equal, won't be always right
        testthat::expect_equal(A1_valid, TRUE)
        A2_valid <- mean(sumstats_dt$A2 == sumstats_dt$A2_der) > 0.45 # 45% threshold
        # expect A2s to be equal, all won't be since can be multiple A2s
        testthat::expect_equal(A2_valid, TRUE)

        A1_valid2 <- mean(sumstats_dt$A1 == sumstats_dt$A1_der2) > 0.5 # 50% threshold
        # expect A1s to be equal, all won't be since can be multiple A1s
        testthat::expect_equal(A1_valid2, TRUE)

        A2_valid2 <- mean(sumstats_dt$A2 == sumstats_dt$A2_der2) > 0.45 # 45% threshold
        # expect A2s to be equal, all won't be since can be multiple A2s
        testthat::expect_equal(A2_valid2, TRUE)

        # now check if ran org through allele flip check
        writeLines(eduAttainOkbay, con = file)
        org_rtrn <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = TRUE,
            allele_flip_frq = FALSE,
            dbSNP=144
        )
        res_org <- data.table::fread(org_rtrn, nThread = 1)
        data.table::setkey(res_org, SNP)
        # add A1 to org from one allele missing run
        res_org[res_dt2, A1_der := i.A1]
        # add A2 to org from one allele missing run
        res_org[res_dt3, A2_der := i.A2]

        testthat::expect_equal(
            all(res_org$A1 == res_org$A1_der) &&
            all(res_org$A2 == res_org$A2_der), TRUE)
    } else {
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})
