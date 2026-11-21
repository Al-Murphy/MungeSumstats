# Regression tests for GitHub issues on reading sumstats (read_vcf, read_sumstats, empty columns).
# Issues: #21, #51, #53, #91/#154, #100, #110/#119, #155/#156.
# read_vcf() writes .bgz/.tbi beside its input, so always read a temp copy
als_tmp <- function(edit = identity) {
    l <- readLines(system.file("extdata", "ALSvcf.vcf",
                               package = "MungeSumstats"))
    f <- file.path(tempfile("als"), "als.vcf")
    dir.create(dirname(f))
    writeLines(edit(l), f)
    f
}
# PGC-style text file: tab-separated, with a VCF-like "#CHROM" header
pgc_tmp <- function(file_name) {
    f <- file.path(tempfile("pgc"), file_name)
    dir.create(dirname(f))
    writeLines(c(
        "#CHROM\tPOS\tID\tA1\tA2\tBETA\tSE\tPVAL\tNGT",
        "1\t10583\trs58108140\tA\tG\t0.01\t0.02\t0.5\t4",
        "1\t30923\trs806731\tG\tT\t-0.02\t0.03\t0.2\t5",
        "2\t51479\trs116400033\tT\tA\t0.03\t0.01\t0.01\t3"
    ), f)
    f
}

test_that("#21: a VCF without ##SAMPLE lines is read (nThread=1)", {
    f <- als_tmp(function(l) l[!startsWith(l, "##SAMPLE")])
    expect_no_error(res <- MungeSumstats::read_vcf(path = f, nThread = 1))
    expect_gt(nrow(res), 0)
    expect_true(all(c("SNP", "ES", "SE", "LP", "AF") %in% names(res)))
})

test_that("#51: a ##SAMPLE ID that differs from the data column is ignored", {
    f <- als_tmp(function(l) {
        sub("^##SAMPLE=<ID=EBI-a-GCST005647", "##SAMPLE=<ID=EBI-a:GCST005647",
            l)
    })
    expect_no_error(res <- MungeSumstats::read_vcf(path = f))
    expect_gt(nrow(res), 0)
    expect_true(all(c("SNP", "ES", "SE", "LP") %in% names(res)))
})

test_that("#100: a VCF INFO 'AF=' field stays AF and is not read as INFO", {
    res <- MungeSumstats::read_vcf(path = als_tmp())
    expect_false("INFO" %in% names(res))
    expect_true("AF" %in% names(res))
    expect_equal(res[SNP == "rs58108140"]$AF[1], 0.1589, tolerance = 1e-6)
})

test_that("#53: read_sumstats accepts an R object", {
    df <- data.frame(MarkerName = c("rs1", "rs2", "rs3"),
                     Pval = c(0.1, 0.2, 0.3))
    expect_message(out <- MungeSumstats::read_sumstats(path = df),
                   "Summary statistics passed as R object")
    expect_true(data.table::is.data.table(out))
    expect_equal(nrow(out), 3)
    expect_equal(nrow(MungeSumstats::read_sumstats(path = df, nrows = 2)), 2)
    # a data.table input is standardised without renaming the caller's object
    dt <- data.table::as.data.table(df)
    std <- MungeSumstats::read_sumstats(path = dt, standardise_headers = TRUE)
    expect_equal(names(std), c("SNP", "P"))
    expect_equal(data.table::key(std), "SNP")
    expect_equal(names(dt), c("MarkerName", "Pval"))
})

test_that("#91/#154: a PGC text file renamed to .tsv is read as tabular", {
    f <- pgc_tmp("pgc.tsv")
    expect_message(out <- MungeSumstats::read_sumstats(path = f),
                   "Importing tabular file")
    # the "#CHROM" header line is kept (it maps to CHR, see test-issues_headers)
    expect_equal(names(out), c("#CHROM", "POS", "ID", "A1", "A2", "BETA",
                               "SE", "PVAL", "NGT"))
    expect_equal(nrow(out), 3)
})

test_that("#91/#154: a PGC '.vcf.tsv' text file is read as tabular", {
    f <- pgc_tmp("pgc.vcf.tsv")
    out <- MungeSumstats::read_sumstats(path = f, standardise_headers = TRUE)
    expect_true(all(c("CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "P",
                      "NGT") %in% names(out)))
    expect_equal(nrow(out), 3)
})

test_that("#155/#156: 'vcf' elsewhere in a tabular path does not route to read_vcf", {
    src <- system.file("extdata", "eduAttainOkbay.txt",
                       package = "MungeSumstats")
    td <- tempfile("vcf_inputs")
    dir.create(file.path(td, "my.vcf.d"), recursive = TRUE)
    paths <- c(file.path(td, "notvcffile.tsv"),
               # "." in ".vcf.tsv" must not act as a regex wildcard
               file.path(td, "gwas_vcf.tsv"),
               file.path(td, "my.vcf.d", "gwas.txt"))
    for (p in paths) {
        file.copy(src, p)
        expect_message(out <- MungeSumstats::read_sumstats(path = p),
                       "Importing tabular file")
        expect_equal(nrow(out), 93)
    }
})

test_that("#110/#119: partly-NA columns are not flagged as empty", {
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"),
        BETA = c(0.1, NA, 0.3),
        INFO = c("0.9", NA, "0.8"),
        `P-value` = c(0.01, NA, 0.5),
        EMPTY_NA = c(NA, NA, NA),
        EMPTY_DOT = c(".", ".", ".")
    )
    ec <- MungeSumstats:::check_empty_cols(sumstats_dt = data.table::copy(dt))
    expect_false(anyNA(names(ec)))
    expect_setequal(names(ec), c("EMPTY_NA", "EMPTY_DOT"))
    # two or more NA names used to give "duplicates detected"
    expect_no_warning(MungeSumstats:::remove_empty_cols(sumstats_dt = dt))
    expect_setequal(names(dt), c("SNP", "BETA", "INFO", "P-value"))
})
