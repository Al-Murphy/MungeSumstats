# Regression tests for GitHub issues on allele flipping (hand-made `rsids`,
# no SNPlocs). Issues: #44, #46, #88, #139, #140, #145, #165.

# rs1: A1 is the reference; rs2: A2 is the reference (flip);
# rs3: neither allele matches; rs4: absent from the reference
flip_ss <- data.table::data.table(
    SNP = c("rs1", "rs2", "rs3", "rs4"), CHR = 1L, BP = 1:4,
    A1 = c("G", "A", "C", "T"), A2 = c("C", "G", "T", "C"),
    BETA = c(0.1, 0.2, 0.3, 0.4), FRQ = c(0.1, 0.2, 0.3, 0.4), P = 0.5
)
flip_rs <- data.table::data.table(
    SNP = c("rs1", "rs2", "rs3"), ref_allele = c("G", "G", "A")
)
data.table::setkey(flip_rs, SNP)

run_flip <- function(ss = flip_ss, rs = flip_rs, ...) {
    args <- utils::modifyList(list(
        sumstats_dt = data.table::copy(ss), path = NULL,
        ref_genome = "GRCh37", rsids = data.table::copy(rs),
        allele_flip_check = TRUE, allele_flip_drop = TRUE,
        allele_flip_z = TRUE, allele_flip_frq = TRUE,
        bi_allelic_filter = TRUE, flip_frq_as_biallelic = FALSE,
        imputation_ind = TRUE, log_folder_ind = FALSE,
        check_save_out = NULL, tabix_index = FALSE, nThread = 1,
        log_files = list(), mapping_file = MungeSumstats:::sumstatsColHeaders,
        dbSNP = 144, dbSNP_tarball = NULL
    ), list(...))
    do.call(MungeSumstats:::check_allele_flip, args)$sumstats_dt
}

test_that("#44/#46/#88: flips invert OR, take 1-FRQ, negate signed stats", {
    # OR-only sumstats (#88: OR used to be negated like BETA)
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2"), CHR = 1L, BP = 1:2,
        A1 = c("G", "A"), A2 = c("C", "G"), OR = c(2, 4), Z = c(1.5, -2),
        LOG_ODDS = c(0.3, 0.4), SIGNED_SUMSTAT = c(1, 2),
        FRQ = c(0.25, 0.1), P = 0.5
    )
    res <- run_flip(ss)
    expect_equal(res[SNP == "rs2", c(A1, A2)], c("G", "A"))
    expect_equal(res[SNP == "rs2", OR], 0.25)
    expect_equal(res[SNP == "rs2", FRQ], 0.9)
    expect_equal(res[SNP == "rs2", Z], 2)
    expect_equal(res[SNP == "rs2", LOG_ODDS], -0.4)
    expect_equal(res[SNP == "rs2", SIGNED_SUMSTAT], -2)
    # flagged when imputation_ind=TRUE; unflipped SNP untouched
    expect_true(res[SNP == "rs2", flipped])
    expect_true(is.na(res[SNP == "rs1", flipped]))
    expect_equal(unlist(res[SNP == "rs1", .(OR, FRQ, Z)]),
                 c(OR = 2, FRQ = 0.25, Z = 1.5))
    # allele_flip_z=FALSE keeps Z; imputation_ind=FALSE adds no flag
    res <- run_flip(ss, allele_flip_z = FALSE, imputation_ind = FALSE)
    expect_equal(res[SNP == "rs2", Z], -2)
    expect_equal(res[SNP == "rs2", OR], 0.25)
    expect_false("flipped" %in% names(res))
})

test_that("#140: allele_flip_drop drops SNPs matching neither allele", {
    expect_message(
        res <- run_flip(),
        paste("There are 1 SNPs where neither A1 nor A2 match the",
              "reference genome. These will be removed.")
    )
    expect_false("rs3" %in% res$SNP)
    expect_equal(res[SNP == "rs2", c(A1, A2)], c("G", "A"))
    expect_equal(res[SNP == "rs2", c(BETA, FRQ)], c(-0.2, 0.8))
    expect_equal(res[SNP == "rs1", c(A1, A2)], c("G", "C"))
    expect_false(any(c("ref_gen_allele", "match_type") %in% names(res)))
    # allele_flip_drop=FALSE keeps them, unflipped; rs2 is still flipped
    expect_no_message(res <- run_flip(allele_flip_drop = FALSE),
                      message = "neither A1 nor A2")
    expect_equal(nrow(res), 4)
    expect_equal(unlist(res[SNP == "rs3", .(A1, A2)]),
                 c(A1 = "C", A2 = "T"))
    expect_equal(res[SNP == "rs3", BETA], 0.3)
    expect_equal(res[SNP == "rs4", c(A1, A2)], c("T", "C"))
    expect_equal(res[SNP == "rs2", BETA], -0.2)
})

test_that("#140: SNPs absent from the reference are kept unflipped", {
    # no SNP matches neither allele, so nothing is dropped
    res <- run_flip(flip_ss[SNP != "rs3"])
    expect_setequal(res$SNP, c("rs1", "rs2", "rs4"))
    expect_equal(res[SNP == "rs4", c(A1, A2)], c("T", "C"))
    expect_equal(res[SNP == "rs4", BETA], 0.4)
    expect_true(is.na(res[SNP == "rs4", flipped]))
})

test_that("#140: allele_flip_drop keeps SNPs absent from the reference", {
    # message reports 1 removal (rs3); rs4 has no reference allele to compare
    res <- run_flip()
    expect_setequal(res$SNP, c("rs1", "rs2", "rs4"))
    expect_equal(res[SNP == "rs4", BETA], 0.4)
})

test_that("#165: FRQ guard names its remedies and only fires on a flip", {
    # the stop message points to both work-arounds
    expect_error(run_flip(bi_allelic_filter = FALSE),
                 "set allele_flip_frq to FALSE")
    expect_error(run_flip(bi_allelic_filter = FALSE),
                 "flip_frq_as_biallelic to TRUE")
    # nothing needs flipping: non-biallelic SNPs kept with FRQ is fine
    expect_no_error(
        res <- run_flip(flip_ss[SNP == "rs1"], bi_allelic_filter = FALSE)
    )
    expect_equal(res$FRQ, 0.1)
    # opt-in reports how many SNPs are flipped as if bi-allelic
    expect_message(
        res <- run_flip(bi_allelic_filter = FALSE,
                        flip_frq_as_biallelic = TRUE),
        "set flip_frq_as_biallelic to TRUE meaning of the 1 SNPs to be flipped"
    )
    expect_equal(res[SNP == "rs2", c(A1, A2)], c("G", "A"))
    expect_equal(res[SNP == "rs2", FRQ], 0.8)
})

test_that("#145: lower-case alleles are uppercased, then match reference", {
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2"), A1 = c("a", "ag"), A2 = c("g", "T")
    )
    expect_message(
        expect_message(
            res <- MungeSumstats:::make_allele_upper(
                sumstats_dt = data.table::copy(ss), log_files = list()
            )$sumstats_dt,
            "Checking A1 is uppercase"
        ),
        "Checking A2 is uppercase"
    )
    expect_equal(res$A1, c("A", "AG"))
    expect_equal(res$A2, c("G", "T"))
    # no allele columns: returned unchanged
    no_al <- data.table::data.table(SNP = "rs1", P = 0.1)
    expect_no_message(
        res <- MungeSumstats:::make_allele_upper(
            sumstats_dt = data.table::copy(no_al), log_files = list()
        )$sumstats_dt
    )
    expect_equal(res, no_al)
    # upper-cased alleles then match the reference instead of all dropping
    lc <- data.table::copy(flip_ss[SNP %in% c("rs1", "rs2")])
    lc[, c("A1", "A2") := list(tolower(A1), tolower(A2))]
    up <- suppressMessages(MungeSumstats:::make_allele_upper(
        sumstats_dt = lc, log_files = list()
    ))$sumstats_dt
    expect_no_message(res <- run_flip(up), message = "neither A1 nor A2")
    expect_equal(nrow(res), 2)
    expect_equal(res[SNP == "rs2", c(A1, A2)], c("G", "A"))
})

test_that("#139: check_bi_allelic keeps off-reference SNPs (none filtered)", {
    # rs3 is an indel absent from the reference; no SNP is non-biallelic
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = 1L, BP = c(10L, 20L, 30L),
        A1 = c("A", "A", "C"), A2 = c("G", "T", "CTT"), BETA = c(0.1, 0.2, 0.3)
    )
    rs <- data.table::data.table(SNP = c("rs1", "rs2"),
                                 alleles_as_ambig = c("R", "W"))
    data.table::setkey(rs, SNP)
    res <- MungeSumstats:::check_bi_allelic(
        sumstats_dt = data.table::copy(ss), path = NULL,
        ref_genome = "GRCh37", bi_allelic_filter = TRUE, rsids = rs,
        log_folder_ind = FALSE, check_save_out = NULL, tabix_index = FALSE,
        nThread = 1, log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
    )$sumstats_dt
    expect_setequal(res$SNP, c("rs1", "rs2", "rs3"))
})

test_that("#139: check_bi_allelic keeps off-reference SNPs (some filtered)", {
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = 1L, BP = c(10L, 20L, 30L),
        A1 = c("A", "A", "C"), A2 = c("G", "T", "CTT"), BETA = c(0.1, 0.2, 0.3)
    )
    # R = biallelic, H = tri-allelic; rs3 (indel) is not in the reference
    rs <- data.table::data.table(SNP = c("rs1", "rs2"),
                                 alleles_as_ambig = c("R", "H"))
    data.table::setkey(rs, SNP)
    expect_message(
        res <- MungeSumstats:::check_bi_allelic(
            sumstats_dt = data.table::copy(ss), path = NULL,
            ref_genome = "GRCh37", bi_allelic_filter = TRUE, rsids = rs,
            log_folder_ind = FALSE, check_save_out = NULL,
            tabix_index = FALSE, nThread = 1, log_files = list(),
            dbSNP = 144, dbSNP_tarball = NULL
        )$sumstats_dt,
        "1 SNPs are non-biallelic"
    )
    expect_setequal(res$SNP, c("rs1", "rs3"))
})
