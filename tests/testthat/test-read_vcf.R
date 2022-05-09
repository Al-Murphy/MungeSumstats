test_that("Test that VCFs can be properly read in.", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        # write the ALS GWAS, VCF file to a temp file for testing
        vcf_path <- system.file("extdata", "ALSvcf.vcf", 
                                package = "MungeSumstats")
        ### Read in original VCF
        sumstats_dt <- read_vcf(path = vcf_path)
        sumstats_dt <- standardise_header(
            sumstats_dt = sumstats_dt,
            return_list = FALSE
        )
        ### Write this VCF and then read it in again
        tmp_vcf <- tempfile(fileext = ".vcf.gz")
        tmp_vcf <- write_sumstats(
            sumstats_dt = sumstats_dt,
            save_path = tmp_vcf,
            write_vcf = TRUE,
            return_path = TRUE, 
            save_path_check = TRUE
        )
        sumstats_dt2 <- read_vcf(path = tmp_vcf)
        sumstats_dt2 <- standardise_header(
            sumstats_dt = sumstats_dt2,
            return_list = FALSE
        )
        ## Reading in the VCF again adds some new cols and changes their order
        common_cols <- intersect(names(sumstats_dt), names(sumstats_dt2))
        testthat::expect_equal(
            sumstats_dt[,common_cols,with=FALSE], 
            sumstats_dt2[,common_cols,with=FALSE]
        )
    }    
    else{
        testthat::expect_equal(is_32bit_windows, TRUE)
    }
})