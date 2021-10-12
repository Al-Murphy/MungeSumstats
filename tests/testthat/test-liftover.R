test_that("liftover", {
    file <- tempfile()
    eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    ))
    # write the Educational Attainment GWAS to a temp file for testing
    writeLines(eduAttainOkbay, con = file)
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" #&&
        #.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        #use stashed versions of chain files to speed up test
        #move them to temp dir so found and avoid download
        save_dir <- tempdir()
        build_conversion <- c("hg38ToHg19", "hg19ToHg38")
        chain_file <- paste0(build_conversion[1],".over.chain.gz")
        chain_file2 <- paste0(build_conversion[2],".over.chain.gz")
        
        
        local_path <- paste0(save_dir,"/",chain_file)
        local_path2 <- paste0(save_dir,"/",chain_file2)
        local_path_gunzip <- gsub(".gz", "", local_path)
        local_path_gunzip2 <- gsub(".gz", "", local_path2)
        
        if(!file.exists(local_path_gunzip)){
            if(!file.exists(local_path))
                copied <- copyFile(srcPathname = 
                                       system.file("extdata",chain_file,
                                                   package="MungeSumstats"), 
                                   save_dir)
            R.utils::gunzip(local_path,overwrite=TRUE)
        }
        if(!file.exists(local_path_gunzip2)){
            if(!file.exists(local_path2))    
                copied2 <- copyFile(srcPathname = 
                                       system.file("extdata",chain_file2,
                                                   package="MungeSumstats"), 
                                   save_dir)
            R.utils::gunzip(local_path2,overwrite=TRUE)
        }
        # Run MungeSumstats code
        reformatted <- MungeSumstats::format_sumstats(
            path = file,
            ref_genome = "GRCh37",
            on_ref_genome = TRUE,
            convert_ref_genome = "GRCh38",
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            log_folder_ind = TRUE,
            imputation_ind = TRUE
        )
        # now rerun and convert back, should then equal original
        reformatted2 <- MungeSumstats::format_sumstats(
            path = reformatted$sumstats,
            ref_genome = "GRCh38",
            on_ref_genome = TRUE,
            convert_ref_genome = "GRCh37",
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            log_folder_ind = TRUE,
            imputation_ind = TRUE
        )

        # run org through
        reformatted3 <- MungeSumstats::format_sumstats(
            path = file,
            ref_genome = "GRCh37",
            on_ref_genome = TRUE,
            convert_ref_genome = NULL,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            log_folder_ind = TRUE
        )

        ref_37_org <- data.table::fread(reformatted3$sumstats, nThread = 1)
        ref_37 <- data.table::fread(reformatted2$sumstats, nThread = 1)
        ref_37[, IMPUTATION_gen_build := NULL]
        ref_37[, IMPUTATION_GEN_BUILD := NULL]
        # drop 1 row from org not in chain files
        ref_37_org <- ref_37_org[SNP %in% ref_37$SNP, ]
        expect_equal(all.equal(ref_37_org, ref_37), TRUE)
    } else {
        expect_equal(is_32bit_windows, TRUE)
    }
})
