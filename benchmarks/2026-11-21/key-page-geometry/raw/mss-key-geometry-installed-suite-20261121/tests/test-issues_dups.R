# Regression tests for GitHub issues on duplicate SNP/BP removal.
# Issues: #77/#81 (bi_allelic_filter/check_dups), #89 and #139 (indels).
run_dup_snp <- function(dt, bi_allelic_filter = TRUE, check_dups = TRUE,
                        indels = TRUE) {
    MungeSumstats:::check_dup_snp(
        sumstats_dt = data.table::copy(dt), indels = indels, path = NULL,
        log_folder_ind = FALSE, check_save_out = list(), tabix_index = FALSE,
        nThread = 1, log_files = list(),
        bi_allelic_filter = bi_allelic_filter, check_dups = check_dups
    )$sumstats_dt
}
run_dup_bp <- function(dt, bi_allelic_filter = TRUE, check_dups = TRUE,
                       indels = TRUE) {
    MungeSumstats:::check_dup_bp(
        sumstats_dt = data.table::copy(dt),
        bi_allelic_filter = bi_allelic_filter, check_dups = check_dups,
        indels = indels, path = NULL, log_folder_ind = FALSE,
        check_save_out = list(), tabix_index = FALSE, nThread = 1,
        log_files = list()
    )$sumstats_dt
}

test_that("#77/#81: RSID dedupe only with bi_allelic_filter and check_dups", {
    # rs1 is tri-allelic: same RSID, different alleles
    dt <- data.table::data.table(
        SNP = c("rs1", "rs1", "rs2"),
        CHR = "1",
        BP = c(100L, 100L, 200L),
        A1 = c("A", "A", "C"),
        A2 = c("G", "T", "T")
    )
    # #81: non-biallelic SNPs kept -> no RSID dedupe
    expect_no_message(
        res <- run_dup_snp(dt, bi_allelic_filter = FALSE),
        message = "RSIDs are duplicated"
    )
    expect_equal(nrow(res), 3)
    expect_message(
        res <- run_dup_snp(dt),
        "1 RSIDs are duplicated in the sumstats file"
    )
    expect_equal(nrow(res), 2)
    # #77: check_dups = FALSE (e.g. QTL sumstats) keeps every row
    expect_equal(nrow(run_dup_snp(dt, check_dups = FALSE)), 3)
})

test_that("#77/#81/#89: BP dedupe skips indels and respects the switches", {
    # rs4 (SNV) and rs34589910 (indel) share 4:6364621
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs4", "rs34589910"),
        CHR = c("1", "1", "4", "4"),
        BP = c(100L, 200L, 6364621L, 6364621L),
        A1 = c("A", "C", "C", "C"),
        A2 = c("G", "T", "G", "CG")
    )
    # #89: the indel is not a duplicate of the SNV at the same position
    expect_message(
        res <- run_dup_bp(dt),
        "Found 1 Indels. These won't be checked for duplicates based on"
    )
    expect_setequal(res$SNP, dt$SNP)
    # indels = FALSE: plain BP dedupe drops the second row at 4:6364621
    expect_message(
        res <- run_dup_bp(dt, indels = FALSE),
        "1 base-pair positions are duplicated in the sumstats file"
    )
    expect_equal(nrow(res), 3)
    # #77/#81: either switch off -> only exact duplicate rows are removed
    res <- suppressMessages(run_dup_bp(dt, check_dups = FALSE))
    expect_equal(nrow(res), 4)
    res <- suppressMessages(run_dup_bp(dt, bi_allelic_filter = FALSE))
    expect_equal(nrow(res), 4)
})

test_that("#139: indels sharing an RSID are deduped on alleles, not RSID", {
    # maintainer's example: one SNV and two identical indels at 8:134625772
    ss <- data.table::data.table(
        SNP = "rs539065912",
        CHR = 8L,
        BP = 134625772L,
        A1 = "G",
        A2 = c("GCCGAAGGAGA", "A", "GCCGAAGGAGA"),
        BETA = c(9.74121e-4, 0.0032263, 6.41437e-4),
        SE = c(0.00753359, 0.0115306, 0.00763662),
        P = c(0.86, 0.99, 0.82)
    )
    expect_message(
        r1 <- run_dup_snp(ss),
        "Found 2 Indels. These won't be checked for duplicates based on RS ID"
    )
    expect_equal(nrow(r1), 3)
    expect_message(
        r2 <- run_dup_bp(r1),
        "1 sumstat rows are duplicated"
    )
    expect_equal(sort(r2$A2), c("A", "GCCGAAGGAGA"))
    # the first indel is kept, as in the maintainer's output
    expect_equal(r2[A2 == "GCCGAAGGAGA"]$BETA, 9.74121e-4)
    # different indel alleles: all three variants kept
    ss2 <- data.table::copy(ss)
    ss2[3, A2 := "GC"]
    res <- suppressMessages(run_dup_bp(run_dup_snp(ss2)))
    expect_equal(nrow(res), 3)
    # indels = FALSE falls back to RSID-only dedupe
    expect_message(
        res <- run_dup_snp(ss, indels = FALSE),
        "2 RSIDs are duplicated"
    )
    expect_equal(nrow(res), 1)
})
