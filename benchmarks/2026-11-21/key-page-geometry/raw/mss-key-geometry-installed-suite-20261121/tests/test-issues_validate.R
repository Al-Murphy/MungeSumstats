# Regression tests for GitHub issues on input-parameter validation.
# Issues: #3, #37, #49, #98, #99, #114, #125, #158, #200, #204.
# Build a full validate_parameters() argument list from the format_sumstats()
# defaults, then apply overrides (a NULL override is kept, not dropped).
vp_args <- function(...) {
    ns <- asNamespace("MungeSumstats")
    fs <- formals(MungeSumstats::format_sumstats)
    fs$path <- system.file("extdata", "eduAttainOkbay.txt",
                           package = "MungeSumstats")
    args <- lapply(fs[names(formals(ns$validate_parameters))], eval,
                   envir = ns)
    utils::modifyList(args, list(...), keep.null = TRUE)
}
validate <- function(...) {
    do.call(MungeSumstats:::validate_parameters, vp_args(...))
}
# validate_parameters() requires the SNPlocs/BSgenome packages of each build.
# Loading a SNPlocs namespace raises harmless "replacing previous import"
# warnings, so they are muffled here.
skip_without <- function(pkgs) {
    for (p in pkgs) suppressWarnings(testthat::skip_if_not_installed(p))
}
refs_37 <- c("SNPlocs.Hsapiens.dbSNP144.GRCh37",
             "BSgenome.Hsapiens.1000genomes.hs37d5")
refs_38 <- c("SNPlocs.Hsapiens.dbSNP144.GRCh38",
             "BSgenome.Hsapiens.NCBI.GRCh38")

test_that("#49/#114: ref_genome=NULL validates; install check uses dbSNP", {
    skip_without(c(refs_37, refs_38,
                   "SNPlocs.Hsapiens.dbSNP155.GRCh37",
                   "SNPlocs.Hsapiens.dbSNP155.GRCh38"))
    # NULL (infer the build) is the format_sumstats() default
    expect_no_error(validate(ref_genome = NULL, dbSNP = 144))
    # #114: with the 155 packages installed, NULL + 155 must not ask for them
    expect_no_error(validate(ref_genome = NULL, dbSNP = 155))
    # the install message names the requested version and build
    expect_error(validate(ref_genome = NULL, dbSNP = 999),
                 "'SNPlocs.Hsapiens.dbSNP999.GRCh37'", fixed = TRUE)
    expect_error(validate(ref_genome = "GRCh38", dbSNP = 999),
                 "'SNPlocs.Hsapiens.dbSNP999.GRCh38'", fixed = TRUE)
})

test_that("#3: ldsc_format is deprecated in favour of save_format", {
    skip_without(refs_37)
    expect_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                          ldsc_format = TRUE),
                 "`ldsc_format` has been deprecated", fixed = TRUE)
    expect_no_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                             save_format = "LDSC"))
})

test_that("#37/#158: rmv_chrPrefix is deprecated; chr_style is validated", {
    skip_without(refs_37)
    expect_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                          rmv_chrPrefix = TRUE),
                 "rmv_chrPrefix is now deprecated")
    expect_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                          chr_style = "bad"),
                 "chr_style must be one of NCBI, UCSC, dbSNP, Ensembl")
    for (s in c("NCBI", "UCSC", "dbSNP", "Ensembl")) {
        expect_no_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                                 chr_style = s))
    }
})

test_that("#98/#200: dbSNP=155 validates; unknown or invalid versions error", {
    skip_without(c("SNPlocs.Hsapiens.dbSNP155.GRCh37",
                   "BSgenome.Hsapiens.1000genomes.hs37d5"))
    expect_no_error(validate(ref_genome = "GRCh37", dbSNP = 155))
    # #200: a missing SNPlocs package gives an actionable install message
    expect_error(
        validate(ref_genome = "GRCh37", dbSNP = 999),
        paste0("To use dbSNP=999 on GRCh37, please install the Bioconductor ",
               "package 'SNPlocs.Hsapiens.dbSNP999.GRCh37'."),
        fixed = TRUE
    )
    expect_error(validate(ref_genome = "GRCh37", dbSNP = -1),
                 "must be a positive integer")
})

test_that("#99: a remote https path passes validation", {
    skip_without(refs_37)
    expect_no_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                             path = "https://example.org/sumstats.txt.gz"))
    # a missing local file is still rejected
    expect_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                          path = tempfile(fileext = ".txt")),
                 "Path to GWAS sumstats is not valid")
})

test_that("#125: several numeric compute_n values give a clear error", {
    skip_without(refs_37)
    expect_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                          compute_n = c(5, 5, 7)),
                 "Multiple integer values for compute_n not supported.",
                 fixed = TRUE)
    # several *methods* are allowed
    expect_no_error(validate(ref_genome = "GRCh37", dbSNP = 144,
                             compute_n = c("ldsc", "giant")))
})

test_that("#204: dbSNP_tarball validation and load_snp_loc_data errors", {
    skip_without("BSgenome.Hsapiens.1000genomes.hs37d5")
    expect_error(
        validate(ref_genome = "GRCh37",
                 dbSNP_tarball = file.path(tempdir(), "nope_GRCh37.tar.gz")),
        "dbSNP_tarball parameter is invalid"
    )
    # an existing tarball replaces the SNPlocs install requirement
    tb <- tempfile("SNPlocs.Hsapiens.dbSNP157.GRCh37_", fileext = ".tar.gz")
    file.create(tb)
    expect_no_error(validate(ref_genome = "GRCh37", dbSNP = 157,
                             dbSNP_tarball = tb))
    unlink(tb)
    # load_snp_loc_data checks the build in the tarball name and the build
    expect_error(
        MungeSumstats::load_snp_loc_data(
            "GRCh38", 156,
            dbSNP_tarball = "SNPlocs.Hsapiens.dbSNP156.GRCh37_0.99.24.tar.gz"
        ),
        "dbSNP_tarball file name does not appear"
    )
    expect_error(
        suppressMessages(MungeSumstats::load_snp_loc_data("GRCh37", 999)),
        "Please install Bioconductor package 'SNPlocs.Hsapiens.dbSNP999.GRCh37'",
        fixed = TRUE
    )
    expect_error(MungeSumstats::load_snp_loc_data("hg19", 155),
                 "must be 'GRCh37' or 'GRCh38'")
})
