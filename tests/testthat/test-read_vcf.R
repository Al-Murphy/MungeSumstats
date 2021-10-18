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
        sumstats_dt <- standardise_sumstats_column_headers_crossplatform(
            sumstats_dt = sumstats_dt,
            mapping_file = sumstatsColHeaders
        )[["sumstats_dt"]]
        ### Write this VCF and then read it in again
        tmp_vcf <- tempfile(fileext = ".vcf.gz")
        write_sumstats(
            sumstats_dt = sumstats_dt,
            sep = "\t",
            save_path = tmp_vcf,
            write_vcf = TRUE
        )
        sumstats_dt2 <- read_vcf(path = tmp_vcf)
        sumstats_dt2 <- standardise_sumstats_column_headers_crossplatform(
            sumstats_dt = sumstats_dt,
            mapping_file = sumstatsColHeaders
        )[["sumstats_dt"]]
        expect_equal(sumstats_dt, sumstats_dt2)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})
