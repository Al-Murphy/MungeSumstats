test_that("liftover", {
    file <- tempfile()
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
                                            package = "MungeSumstats"
    ))
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay, con = file)
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" &&
        .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(
            path = file,
              ref_genome = "GRCh37",
              on_ref_genome = TRUE,
              convert_ref_genome="GRCh38",
              strand_ambig_filter = FALSE,
              bi_allelic_filter = FALSE,
              allele_flip_check = FALSE,
              log_folder_ind = TRUE,
              imputation_ind = TRUE)
        # now rerun and convert back, should then equal original
        reformatted2 <- MungeSumstats::format_sumstats(
            path = reformatted$sumstats,
              ref_genome = "GRCh38",
              on_ref_genome = TRUE,
              convert_ref_genome="GRCh37",
              strand_ambig_filter = FALSE,
              bi_allelic_filter = FALSE,
              allele_flip_check = FALSE,
              log_folder_ind = TRUE,
              imputation_ind = TRUE)
        
        #run org through
        reformatted3 <- MungeSumstats::format_sumstats(
            path = file,
              ref_genome = "GRCh37",
              on_ref_genome = TRUE,
              convert_ref_genome=NULL,
              strand_ambig_filter = FALSE,
              bi_allelic_filter = FALSE,
              allele_flip_check = FALSE,
              log_folder_ind = TRUE)
        
        ref_37_org <- data.table::fread(reformatted3$sumstats, nThread = 1)
        ref_37 <- data.table::fread(reformatted2$sumstats, nThread = 1)
        ref_37[,IMPUTATION_gen_build:=NULL]
        ref_37[,IMPUTATION_GEN_BUILD:=NULL]
        #drop 1 row from org not in chain files
        ref_37_org <- ref_37_org[SNP %in% ref_37$SNP,]
        expect_equal(all.equal(ref_37_org,ref_37),TRUE)
        
    } else {
        expect_equal(is_32bit_windows, TRUE)
    }
})
