# Long-running regression tests for GitHub issues that need full SNPlocs/
# BSgenome lookups. Issues: #28, #49, #98, #153.
# Genome-build inference looks up both builds (dbSNP 144).
skip_without_both_builds <- function() {
    pkgs <- c("SNPlocs.Hsapiens.dbSNP144.GRCh37",
              "SNPlocs.Hsapiens.dbSNP144.GRCh38",
              "BSgenome.Hsapiens.1000genomes.hs37d5",
              "BSgenome.Hsapiens.NCBI.GRCh38")
    for (p in pkgs) suppressWarnings(testthat::skip_if_not_installed(p))
}

test_that("#49: format_sumstats(ref_genome=NULL) infers GRCh37", {
    skip_without_both_builds()
    file <- tempfile(fileext = ".txt")
    file.copy(system.file("extdata", "eduAttainOkbay.txt",
                          package = "MungeSumstats"), file)
    expect_message(
        res <- MungeSumstats::format_sumstats(
            path = file, ref_genome = NULL, dbSNP = 144,
            on_ref_genome = FALSE, strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE, allele_flip_check = FALSE,
            return_data = TRUE
        ),
        "Inferred genome build: GRCH37"
    )
    expect_equal(nrow(res), 93)
    expect_true(all(c("SNP", "CHR", "BP", "A1", "A2", "P") %in% names(res)))
})

test_that("#153: get_genome_build strips 'chr' and errors on no matches", {
    skip_without_both_builds()
    fe <- MungeSumstats::formatted_example()
    chr_dt <- data.table::copy(fe)[, CHR := paste0("chr", CHR)]
    expect_equal(
        suppressMessages(MungeSumstats:::get_genome_build(
            sumstats = chr_dt, standardise_headers = FALSE, dbSNP = 144
        )),
        "GRCH37"
    )
    # nothing matches either build: error instead of defaulting to GRCh38
    bad <- data.table::copy(fe)[, BP := 1L]
    expect_error(
        suppressMessages(MungeSumstats:::get_genome_build(
            sumstats = bad, standardise_headers = FALSE, dbSNP = 144
        )),
        "No matches found in either reference genome"
    )
})

test_that("#28: get_genome_build ignores NA and '.' SNP IDs", {
    skip_without_both_builds()
    fe <- MungeSumstats::formatted_example()
    fe[1:5, SNP := NA]
    fe[6:8, SNP := "."]
    expect_equal(
        suppressMessages(MungeSumstats:::get_genome_build(
            sumstats = fe, standardise_headers = FALSE, dbSNP = 144
        )),
        "GRCH37"
    )
})

test_that("#98: dbSNP 155 SNPlocs data loads and resolves an rsID", {
    suppressWarnings(
        testthat::skip_if_not_installed("SNPlocs.Hsapiens.dbSNP155.GRCh37")
    )
    expect_message(
        snp_locs <- MungeSumstats::load_snp_loc_data(ref_genome = "GRCh37",
                                                     dbSNP = 155),
        "Loading SNPlocs data for build 155 on GRCH37"
    )
    expect_s4_class(snp_locs, "ODLT_SNPlocs")
    gpos <- BSgenome::snpsById(snp_locs, "rs301800", ifnotfound = "drop")
    expect_equal(as.character(GenomicRanges::seqnames(gpos)), "1")
    expect_equal(GenomicRanges::start(gpos), 8490603L)
})
