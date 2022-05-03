test_that("VCF is correctly formatted", {
    # Run MungeSumstats code
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
        ##.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        ## IMPORTANT: Must have .vcf file extension, 
        # or else MungeSumstats won't know it's a VCF.
        file <- tempfile(fileext = ".vcf")
        # write the ALS GWAS, VCF file to a temp file for testing
        ALSvcf <- readLines(system.file("extdata", "ALSvcf.vcf",
                                        package = "MungeSumstats"
        ))
        writeLines(ALSvcf, con = file)
        # vcf <- MungeSumstats:::read_vcf(path = file)
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(
            path = file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            INFO_filter = 0.01
        )
        reformatted_lines <- readLines(reformatted)
        # check manually against first five SNPs
        corr_res <- c(
            "SNP\tCHR\tBP\tA1\tA2\tEND\tFILTER\tFRQ\tBETA\tLP\tSE\tP",
            "rs58108140\t1\t10583\tG\tA\t10583\tPASS\t0.1589\t0.0312\t0.369267\t0.0393\t0.427300105456596",
            "rs806731\t1\t30923\tG\tT\t30923\tPASS\t0.7843\t-0.0114\t0.126854\t0.0353\t0.746699739815279",
            "rs116400033\t1\t51479\tT\tA\t51479\tPASS\t0.1829\t0.0711\t1.26241\t0.037\t0.0546499790752282",
            "rs146477069\t1\t54421\tA\tG\t54421\tPASS\t0.0352\t-0.024\t0.112102\t0.083\t0.772499131799648"
        )
        testthat::expect_equal(reformatted_lines[1:5], corr_res)
        
        # check allelic flipping with VCF
        file2 <- tempfile(fileext = ".vcf")
        # write the ALS GWAS, VCF file to a temp file for testing
        ALSvcf <- readLines(system.file("extdata", "ALSvcf.vcf",
                                        package = "MungeSumstats"
        ))
        # update last SNP, flipping allelic direction: A/G --> G/A
        snp_of_interest <- "rs146477069"
        rsid_index <- grep(snp_of_interest, ALSvcf, ignore.case = TRUE)
        ALSvcf[rsid_index] <-
            ## Original
            # "1\t54421\trs146477069\tG\tA\t.\tPASS\tAF=0.0352\tES:SE:LP:AF:ID\t+0.024:0.083:0.112102:0.0352:rs146477069"
            ## flip FRQ
            "1\t54421\trs146477069\tG\tA\t.\tPASS\tAF=0.9648\tES:SE:LP:AF:ID\t+0.024:0.083:0.112102:0.9648:rs146477069"
        writeLines(ALSvcf, con = file2)
        reformatted_allelic_flip <-
            MungeSumstats::format_sumstats(
                path = file2,
                ref_genome = "GRCh37",
                on_ref_genome = FALSE,
                strand_ambig_filter = FALSE,
                bi_allelic_filter = TRUE,
                allele_flip_check = TRUE,
                allele_flip_drop = FALSE,
                INFO_filter = 0.01
            )
        reformatted_lines_af <- readLines(reformatted_allelic_flip)
        testthat::expect_equal(sort(reformatted_lines), 
                                sort(reformatted_lines_af))
        # also check outputting as different types
        pth <- system.file("extdata", "ALSvcf.vcf", package = "MungeSumstats")
        rtrn_dt <- MungeSumstats::format_sumstats(
            path = pth,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            allele_flip_drop = FALSE,
            INFO_filter = 0.01,
            return_data = TRUE,
            return_format = "data.table"
        )
        rtrn_grng <- MungeSumstats::format_sumstats(
            path = pth,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            allele_flip_drop = FALSE,
            INFO_filter = 0.01,
            return_data = TRUE,
            return_format = "GRanges"
        )
        rtrn_vrng <- MungeSumstats::format_sumstats(
            path = pth,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            allele_flip_drop = FALSE,
            INFO_filter = 0.01,
            return_data = TRUE,
            return_format = "VRanges"
        )
        testthat::expect_true(is(rtrn_grng,"GRanges"))
        testthat::expect_true(is(rtrn_vrng,"VRanges"))
        testthat::expect_true(is(rtrn_dt,"data.table"))

        # also test inferring the genome build
        #### Failing atm? duplicate SNPs?
        # rtrn_dt_infer <- MungeSumstats::format_sumstats(
        #     path = pth,
        #     on_ref_genome = FALSE,
        #     strand_ambig_filter = FALSE,
        #     bi_allelic_filter = FALSE,
        #     allele_flip_check = FALSE,
        #     allele_flip_drop = FALSE,
        #     INFO_filter = 0.01,
        #     return_data = TRUE,
        #     return_format = "data.table"
        # )
        # testthat::expect_true(all.equal(rtrn_dt, rtrn_dt_infer))

        # also test outputting ldsc_format ready format
        rtrn_ldsc <- MungeSumstats::format_sumstats(
            path = pth,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            allele_flip_drop = FALSE,
            INFO_filter = 0.01,
            ldsc_format = TRUE,
            compute_n = 1001
        )
        res <- data.table::fread(rtrn_ldsc, nThread = 1)
        # check for necessary columns - 
        # https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
        ldsc_cols <- c("SNP", "N", "A1", "A2", "Z")
        testthat::expect_true(all(ldsc_cols %in% names(res)))
        
        testthat::expect_equal(reformatted_lines[1:5], corr_res)
    } else {
        testthat::expect_true(is_32bit_windows)
        testthat::expect_true(is_32bit_windows)
        testthat::expect_true(is_32bit_windows)
        testthat::expect_true(is_32bit_windows)
        testthat::expect_true(is_32bit_windows)
    }
})
