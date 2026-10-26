test_that("Handle missing data", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # Remove data from line 3 to check it is deleted
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
            package = "MungeSumstats"
        ))
        eduAttainOkbay_missing <- eduAttainOkbay
        eduAttainOkbay_missing[3] <-
            "rs12987662\t2\t100821548\tA\tC\t0.3787\t0.027\t0.003\t"
        problem_snp <- "rs9320913"
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay_missing, con = file)
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            sort_coordinates = FALSE,
            dbSNP=144
        )
        reformatted_lines <- readLines(reformatted)
        # Should equal org apart from this one line
        writeLines(eduAttainOkbay, con = file)
        org <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            on_ref_genome = FALSE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            sort_coordinates = FALSE,
            dbSNP=144
        )
        org_lines <- readLines(org)
        rsid_index <- grep(problem_snp, org_lines, ignore.case = TRUE)
        # reordering in function, line 3 is now 58
        expect_equal(reformatted_lines, org_lines[-rsid_index])
        # test imputing CHR BP from SNP when they are na but cols exist
        miss <- fread(system.file("extdata", "eduAttainOkbay.txt",
                                  package = "MungeSumstats"
        ))
        #add NA's
        miss[MarkerName=='rs9320913',CHR:=NA]
        miss[MarkerName=='rs9320913',POS:=NA]
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(miss,
                                            ref_genome = "GRCh37",
                                            on_ref_genome = FALSE,
                                            strand_ambig_filter = FALSE,
                                            bi_allelic_filter = FALSE,
                                            allele_flip_check = FALSE,
                                            sort_coordinates = FALSE,
                                            dbSNP=144
        )
        reformatted_lines <- readLines(reformatted)
        testthat::expect_equal(reformatted_lines, org_lines)
        
        # set `drop_na_cols` to `NULL`
        miss_extra_col <- miss
        miss_extra_col$extra <- NA
        
        testthat::expect_error(MungeSumstats::format_sumstats(
          miss_extra_col,
          ref_genome = "GRCh37",
          on_ref_genome = FALSE,
          strand_ambig_filter = FALSE,
          bi_allelic_filter = FALSE,
          allele_flip_check = FALSE,
          sort_coordinates = FALSE,
          dbSNP = 144, 
          drop_na_cols = NULL
        ), 
        regexp = "All SNPs have been filtered out of  your summary statistics dataset")
          
        reformatted_extra_col <- MungeSumstats::format_sumstats(
          miss_extra_col,
          ref_genome = "GRCh37",
          on_ref_genome = FALSE,
          strand_ambig_filter = FALSE,
          bi_allelic_filter = FALSE,
          allele_flip_check = FALSE,
          sort_coordinates = FALSE,
          dbSNP = 144, 
          drop_na_cols = c("CHRA", "APOS")
        )
        
        reformatted_extra_col_lines <- readLines(reformatted_extra_col)
        expect_equal(length(reformatted_extra_col_lines), length(org_lines))
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})
