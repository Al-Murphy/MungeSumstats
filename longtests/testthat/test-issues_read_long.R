# Slow regression tests for GitHub issues on reading VCFs with several threads
# (spawns BiocParallel workers). Issues: #21.
# ALSvcf.vcf header without ##SAMPLE (so TotalVariants is unknown), written to
# a fresh dir because read_vcf() writes .bgz/.tbi beside its input
als_no_sample <- function(all_autosomes) {
    l <- readLines(system.file("extdata", "ALSvcf.vcf",
                               package = "MungeSumstats"))
    hdr <- l[startsWith(l, "#") & !startsWith(l, "##SAMPLE")]
    body <- unique(l[!startsWith(l, "#")])
    if (all_autosomes) {
        # one record per autosome (all records are on chr1 in the original)
        body <- vapply(seq_len(22), function(i) {
            sub("^1\t", paste0(i, "\t"), body[i])
        }, character(1))
    }
    f <- file.path(tempfile("als"), "als.vcf")
    dir.create(dirname(f))
    writeLines(c(hdr, body), f)
    f
}

test_that("#21: read_vcf(nThread=2) reads a VCF without ##SAMPLE lines", {
    skip_if_not_installed("GenomicFiles")
    skip_if_not_installed("BiocParallel")
    expect_no_error(
        mt <- MungeSumstats::read_vcf(path = als_no_sample(TRUE), nThread = 2)
    )
    st <- MungeSumstats::read_vcf(path = als_no_sample(TRUE), nThread = 1)
    expect_equal(nrow(mt), 22)
    data.table::setkeyv(mt, "SNP")
    data.table::setkeyv(st, "SNP")
    expect_equal(mt[, names(st), with = FALSE], st, ignore_attr = TRUE)
})

test_that("#21: read_vcf(nThread=2) skips header contigs that have no records", {
    skip_if_not_installed("GenomicFiles")
    skip_if_not_installed("BiocParallel")
    expect_no_error(
        res <- MungeSumstats::read_vcf(path = als_no_sample(FALSE), nThread = 2)
    )
    expect_equal(nrow(res), 101)
})
