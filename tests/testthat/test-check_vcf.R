test_that("Checking whether a file is works", {
    file <- tempfile()
    # write the ALS GWAS, VCF file to a temp file for testing
    vcf_head <- read_header(path = system.file("extdata", "ALSvcf.vcf",
        package = "MungeSumstats"
    ))
    is_vcf <- check_vcf(header = vcf_head)
    expect_equal(is_vcf, TRUE)
})
