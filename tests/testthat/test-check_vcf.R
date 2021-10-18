test_that("Checking whether a file is works", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        # write the ALS GWAS, VCF file to a temp file for testing
        vcf_head <- read_header(path = system.file("extdata", "ALSvcf.vcf",
            package = "MungeSumstats"
        ))
        is_vcf <- check_vcf(header = vcf_head)
        expect_equal(is_vcf, TRUE)
    }    
    else{
        expect_equal(TRUE, TRUE)
    }
})
