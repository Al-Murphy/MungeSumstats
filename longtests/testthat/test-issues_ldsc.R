# Regression tests for GitHub issues on LDSC output allele naming (MSS: A2 is
# the effect allele; LDSC: A1 is). Issues: #152, #160.

test_that("#152/#160: save_format='LDSC' renames A1/A2, effects not negated", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    skip_if_not_installed("BSgenome.Hsapiens.1000genomes.hs37d5")
    ss <- MungeSumstats::formatted_example()[1:20]
    ss[, N := 1000L]
    args <- list(
        ref_genome = "GRCh37", dbSNP = 144, on_ref_genome = FALSE,
        strand_ambig_filter = FALSE, bi_allelic_filter = TRUE,
        allele_flip_check = TRUE, compute_z = "BETA", nThread = 1,
        return_data = TRUE
    )
    mss <- do.call(MungeSumstats::format_sumstats,
                   c(list(path = data.table::copy(ss)), args))
    expect_message(
        ldsc <- do.call(MungeSumstats::format_sumstats,
                        c(list(path = data.table::copy(ss),
                               save_format = "LDSC"), args)),
        "Renaming A1,A2 to match LDSC format."
    )
    data.table::setkey(mss, SNP)
    data.table::setkey(ldsc, SNP)
    orig <- ss[match(mss$SNP, ss$SNP)]
    # MSS output: A1 = reference, A2 = effect; SNPs whose alleles were
    # swapped to get there (some exist in this input) had BETA negated
    swapped <- mss$A1 != orig$A1
    expect_true(any(swapped))
    expect_equal(mss$A1[swapped], orig$A2[swapped])
    expect_equal(mss$BETA[swapped], -orig$BETA[swapped])
    expect_equal(mss$BETA[!swapped], orig$BETA[!swapped])
    # LDSC output: same SNPs, A1/A2 swapped by name, effects unchanged
    expect_equal(ldsc$SNP, mss$SNP)
    expect_equal(ldsc$A1, mss$A2)
    expect_equal(ldsc$A2, mss$A1)
    expect_equal(ldsc$BETA, mss$BETA)
    expect_equal(ldsc$Z, mss$Z)
    expect_equal(ldsc$FRQ, mss$FRQ)
})
