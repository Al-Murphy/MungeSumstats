test_that("strand-ambiguous SNPs are removed", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
        ##.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # Remove all known strand ambiguous SNPs from dataset
        # "rs7131944"  "rs12682297" "rs55830725" "rs2992632"  "rs11689269"
        # "rs4493682"  "rs34106693" "rs6799130"  are strand-ambiguous
        SA_snps <- c(
            "rs7131944", "rs12682297", "rs55830725", "rs2992632", "rs11689269",
            "rs4493682", "rs34106693", "rs6799130"
        )
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
                                                package = "MungeSumstats"
        ))
        eduAttainOkbay_missing <- eduAttainOkbay
        for (SA_snp in SA_snps) {
            eduAttainOkbay_missing <-
                eduAttainOkbay_missing[!grepl(SA_snp, eduAttainOkbay_missing)]
        }
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay_missing, con = file)
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = TRUE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        reformatted_lines <- readLines(reformatted)
        # Should equal org as strand ambig can be removed
        writeLines(eduAttainOkbay, con = file)
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = TRUE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            dbSNP=144
        )
        org_lines <- readLines(org)
        # reordering in function, line 3 rs9320913 is now 58
        expect_equal(setequal(reformatted_lines, org_lines), TRUE)
    } else {
        expect_equal(is_32bit_windows, TRUE)
    }
})
