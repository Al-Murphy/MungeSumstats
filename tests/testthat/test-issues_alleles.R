# Regression tests for GitHub issues on allele handling.
# Hand-made `rsids` reference tables are passed so no SNPlocs/BSgenome is
# loaded. Issues: #111/#124 (bi-allelic filter), #137 (allele flipping).
test_that("#111/#124: bi_allelic_filter drops >2 allele SNPs, FALSE keeps", {
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = "1", BP = 1:3,
        A1 = c("A", "A", "C"), A2 = c("G", "C", "T"), P = 0.1
    )
    # IUPAC: R = A/G, V = A/C/G (tri-allelic), Y = C/T
    rs <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"),
        alleles_as_ambig = c("R", "V", "Y")
    )
    data.table::setkey(rs, SNP)
    run_bi <- function(filt, rsids) {
        MungeSumstats:::check_bi_allelic(
            sumstats_dt = data.table::copy(dt), path = NULL,
            ref_genome = "GRCh37", bi_allelic_filter = filt, rsids = rsids,
            log_folder_ind = FALSE, check_save_out = NULL,
            tabix_index = FALSE, nThread = 1, log_files = list(),
            dbSNP = 144, dbSNP_tarball = NULL
        )
    }
    # filter off: nothing removed and no reference loaded
    expect_no_message(res <- run_bi(FALSE, NULL),
                      message = "Checking for bi-allelic SNPs")
    expect_equal(nrow(res$sumstats_dt), 3)
    expect_null(res$rsids)
    # filter on: tri-allelic rs2 removed
    expect_message(res <- run_bi(TRUE, rs),
                   "1 SNPs are non-biallelic. These will be removed.")
    expect_setequal(res$sumstats_dt$SNP, c("rs1", "rs3"))
})

test_that("#137: allele flip contract and non-biallelic FRQ guard", {
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = "1", BP = 1:3,
        A1 = c("A", "A", "C"), A2 = c("G", "G", "G"),
        BETA = c(0.1, 0.2, 0.3), Z = c(1, 2, 3), OR = c(1.1, 2, 1.2),
        FRQ = c(0.2, 0.3, 0.4), P = 0.01
    )
    # rs1 matches, rs2 needs flipping, rs3 matches neither allele
    rs <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"),
        ref_allele = c("A", "G", "T")
    )
    data.table::setkey(rs, SNP)
    run_flip <- function(flip_frq, bi_allelic, frq_as_bi) {
        MungeSumstats:::check_allele_flip(
            sumstats_dt = data.table::copy(dt), path = NULL,
            ref_genome = "GRCh37", rsids = data.table::copy(rs),
            allele_flip_check = TRUE, allele_flip_drop = TRUE,
            allele_flip_z = TRUE, allele_flip_frq = flip_frq,
            bi_allelic_filter = bi_allelic,
            flip_frq_as_biallelic = frq_as_bi,
            imputation_ind = FALSE, log_folder_ind = FALSE,
            check_save_out = NULL, tabix_index = FALSE, nThread = 1,
            log_files = list(),
            mapping_file = MungeSumstats:::sumstatsColHeaders,
            dbSNP = 144, dbSNP_tarball = NULL
        )$sumstats_dt
    }
    # flipping FRQ while keeping non-biallelic SNPs is refused
    expect_error(run_flip(TRUE, FALSE, FALSE), "only bi-allelic SNPs can be")
    # no FRQ flip: alleles and effects flipped, FRQ untouched
    expect_message(res <- run_flip(FALSE, FALSE, FALSE),
                   "There are 1 SNPs where A1 doesn't match")
    expect_equal(nrow(res), 2)
    expect_equal(res[SNP == "rs2", c(A1, A2)], c("G", "A"))
    expect_equal(res[SNP == "rs2", BETA], -0.2)
    expect_equal(res[SNP == "rs2", Z], -2)
    expect_equal(res[SNP == "rs2", OR], 0.5)
    expect_equal(res[SNP == "rs2", FRQ], 0.3)
    expect_equal(res[SNP == "rs1", c(A1, A2)], c("A", "G"))
    expect_equal(res[SNP == "rs1", BETA], 0.1)
    # bi-allelic data: FRQ is flipped to 1-FRQ
    res <- run_flip(TRUE, TRUE, FALSE)
    expect_equal(res[SNP == "rs2", FRQ], 0.7)
    expect_equal(res[SNP == "rs2", BETA], -0.2)
    # user asserts the SNPs are bi-allelic
    expect_message(res <- run_flip(TRUE, FALSE, TRUE),
                   "flip_frq_as_biallelic")
    expect_equal(res[SNP == "rs2", FRQ], 0.7)
})
