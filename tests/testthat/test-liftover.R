test_that("liftover works", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" #&&
        #.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        ##### Test exported function directly ####
        sumstats_dt <- MungeSumstats::formatted_example()
        
        ## as data.table  
        sumstats_dt_hg38 <- MungeSumstats::liftover(
            sumstats_dt=sumstats_dt,
            ref_genome = "hg19",
            convert_ref_genome="hg38")
        testthat::expect_true("IMPUTATION_gen_build" %in% colnames(sumstats_dt_hg38))
        testthat::expect_equal(nrow(sumstats_dt_hg38), nrow(sumstats_dt))
        proportion_positions_changed <- sum(sumstats_dt_hg38$BP!=sumstats_dt$BP)/
            nrow(sumstats_dt)
        testthat::expect_equal(proportion_positions_changed,1)   
        proportion_snps_changed <- sum(sumstats_dt_hg38$SNP!=sumstats_dt$SNP)/
            nrow(sumstats_dt)
        testthat::expect_equal(proportion_snps_changed,0)     
        
        ## as GRAnges  
        granges_hg38 <- MungeSumstats::liftover(
            sumstats_dt=sumstats_dt,
            ref_genome = "hg19",
            convert_ref_genome="hg38", 
            as_granges = TRUE)
        testthat::expect_true("IMPUTATION_gen_build" %in% colnames(sumstats_dt_hg38))
        testthat::expect_equal(nrow(sumstats_dt_hg38), nrow(sumstats_dt))
        proportion_positions_changed <- sum(GenomicRanges::start(granges_hg38)!=sumstats_dt$BP)/
            nrow(sumstats_dt)
        testthat::expect_equal(proportion_positions_changed,1)   
        proportion_snps_changed <- sum(granges_hg38$SNP!=sumstats_dt$SNP)/
            nrow(sumstats_dt)
        testthat::expect_equal(proportion_snps_changed,0)     
        
        
        
        
        #### Test liftover functionality within format_sumstats ####
        
        file <- tempfile()
        eduAttainOkbay <- readLines(system.file("extdata", "eduAttainOkbay.txt",
                                                package = "MungeSumstats"
        ))
        # write the Educational Attainment GWAS to a temp file for testing
        writeLines(eduAttainOkbay, con = file)
        #use stashed versions of chain files to speed up test
        #move them to temp dir so found and avoid download
        save_dir <- tempdir()
        build_conversion <- c("GRCh38_to_GRCh37", "GRCh37_to_GRCh38")
        chain_file <- paste0(build_conversion[1],".chain.gz")
        chain_file2 <- paste0(build_conversion[2],".chain.gz")
        
        
        local_path <- paste0(save_dir,"/",chain_file)
        local_path2 <- paste0(save_dir,"/",chain_file2)
        local_path_gunzip <- gsub(".gz", "", local_path)
        local_path_gunzip2 <- gsub(".gz", "", local_path2)
        
        if(!file.exists(local_path_gunzip)){
            if(!file.exists(local_path))
                copied <- R.utils::copyFile(srcPathname = 
                                       system.file("extdata",chain_file,
                                                   package="MungeSumstats"), 
                                   save_dir)
            R.utils::gunzip(local_path,overwrite=TRUE)
        }
        if(!file.exists(local_path_gunzip2)){
            if(!file.exists(local_path2))    
                copied2 <- R.utils::copyFile(srcPathname = 
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
            imputation_ind = TRUE,
            dbSNP=144
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
            imputation_ind = TRUE,
            dbSNP=144
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
            log_folder_ind = TRUE,
            dbSNP=144
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
