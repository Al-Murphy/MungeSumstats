# Long regression tests for GitHub issues that need full reference lookups.
# Issues: #23 (junk IDs before snpsById), #83 (chr:bp IDs with NA CHR/BP).
test_that("#23: load_ref_genome_data drops NA, '.' and multi-rs IDs", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    skip_if_not_installed("BSgenome.Hsapiens.1000genomes.hs37d5")
    snps <- c("rs12987662", "rs9320913", NA, ".", "rs5772025_rs397784053")
    # these used to reach BSgenome::snpsById -> "Error in ids2rowids(ids)"
    # (junk IDs raise "NAs introduced by coercion" warnings by design)
    expect_no_error(out <- suppressWarnings(suppressMessages(
        MungeSumstats:::load_ref_genome_data(
            snps = snps,
            ref_genome = "GRCh37",
            dbSNP = 144
        )
    )))
    expect_equal(nrow(out), 2)
    expect_setequal(out$SNP, c("rs12987662", "rs9320913"))
})

test_that("#83: chr:bp_allele SNP IDs with NA CHR/BP don't crash IRanges", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    skip_if_not_installed("BSgenome.Hsapiens.1000genomes.hs37d5")
    # maintainer's reproducer from the issue thread
    a <- data.table::fread(system.file("extdata", "eduAttainOkbay.txt",
                                       package = "MungeSumstats"))
    a[MarkerName %in% c("rs12987662", "rs9320913"), POS := NA]
    a[MarkerName %in% c("rs12987662", "rs9320913"), CHR := NA]
    a[MarkerName == "rs12987662", MarkerName := "2:100821548_AC_A"]
    a[MarkerName == "rs9320913", MarkerName := "6:98584733_AC_A"]
    expect_no_error(out <- suppressWarnings(suppressMessages(
        MungeSumstats::format_sumstats(
            path = a,
            ref_genome = "GRCh37",
            dbSNP = 144,
            on_ref_genome = TRUE,
            strand_ambig_filter = FALSE,
            bi_allelic_filter = FALSE,
            allele_flip_check = FALSE,
            return_data = TRUE
        )
    )))
    # each problem ID is imputed or dropped, never kept as chr:bp
    expect_true(all(grepl("^rs", out$SNP)))
    expect_false(anyNA(out$CHR))
    expect_false(anyNA(out$BP))
    expect_gte(nrow(out), 89)
})
