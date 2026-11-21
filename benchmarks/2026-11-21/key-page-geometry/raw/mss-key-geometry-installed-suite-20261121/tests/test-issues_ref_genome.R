# Regression tests for GitHub issues on reference-genome checks/imputation.
# Issues: #164, #166 (check_on_ref_genome), #150 (NA CHR/BP), #139, #177.
keyed_rsids <- function(...) {
    rs <- data.table::data.table(...)
    data.table::setkey(rs, SNP)
    rs
}
run_on_ref <- function(dt, rsids, indels = TRUE, imputation_ind = FALSE) {
    MungeSumstats:::check_on_ref_genome(
        sumstats_dt = data.table::copy(dt), path = "",
        ref_genome = "GRCh37", on_ref_genome = TRUE, indels = indels,
        rsids = data.table::copy(rsids), imputation_ind = imputation_ind,
        log_folder_ind = FALSE, check_save_out = list(), tabix_index = FALSE,
        nThread = 1, log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
    )
}

# rs2000000000 is not in dbSNP144; its GRCh37 CHR:BP is rs11712056
off_ref_dt <- data.table::data.table(
    SNP = c("rs12987662", "rs9320913", "rs2000000000"),
    CHR = c(2L, 6L, 3L),
    BP = c(100821548L, 98584733L, 49914397L),
    A1 = c("A", "A", "T"),
    A2 = c("C", "C", "C"),
    BETA = c(0.027, 0.024, 0.024),
    SE = 0.003,
    P = c(2.693e-24, 2.457e-19, 3.304e-19)
)
on_ref <- keyed_rsids(SNP = c("rs12987662", "rs9320913"))

test_that("#164: many indels sharing one ID don't explode the not-join", {
    # off-reference indels sharing an ID: the not-join on a duplicated
    # key built n x n row indices (2.5e9 > 2^31 for n = 50,000)
    no_ref <- keyed_rsids(SNP = character(0))
    n <- 50000L
    dt <- data.table::data.table(
        SNP = rep("rs999", n), CHR = 1L, BP = 1000L,
        A1 = "A", A2 = "AT", BETA = 0.1, P = 0.5
    )
    expect_message(res <- run_on_ref(dt, no_ref), "Found 50,000 Indels")
    expect_equal(nrow(res$sumstats_dt), n)
    # small companion: off-reference indels are kept, not dropped
    expect_message(res <- run_on_ref(dt[1:3], no_ref), "Found 3 Indels")
    expect_equal(nrow(res$sumstats_dt), 3)
})

test_that("#166: off-reference SNP with imputation_ind=TRUE, indels T/F", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    # skip the slow snpsById reload; check_no_snp's position lookup is real
    testthat::local_mocked_bindings(
        load_ref_genome_data = function(...) data.table::copy(on_ref),
        .package = "MungeSumstats"
    )
    for (indels in c(FALSE, TRUE)) {
        # indels=FALSE used to rbind a 0-column indel_dt -> column mismatch
        expect_message(
            res <- run_on_ref(off_ref_dt, on_ref, indels = indels,
                              imputation_ind = TRUE)$sumstats_dt,
            "1 SNPs are not on the reference genome"
        )
        expect_equal(nrow(res), 3)
        expect_setequal(res$SNP,
                        c("rs12987662", "rs9320913", "rs11712056"))
        expect_true("IMPUTATION_SNP" %in% names(res))
        expect_equal(res[IMPUTATION_SNP %in% TRUE]$SNP, "rs11712056")
    }
})

test_that("#164: off-reference SNVs are still corrected when indels exist", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    testthat::local_mocked_bindings(
        load_ref_genome_data = function(...) data.table::copy(on_ref),
        .package = "MungeSumstats"
    )
    run_imp <- function(dt) {
        run_on_ref(dt, on_ref, indels = TRUE,
                   imputation_ind = TRUE)$sumstats_dt
    }
    # a) an unrelated off-reference indel hides the off-reference SNV
    a <- rbind(off_ref_dt, data.table::data.table(
        SNP = "rs3", CHR = 4L, BP = 100L, A1 = "A", A2 = "AT",
        BETA = 0.1, SE = 0.003, P = 0.5
    ))
    expect_message(res <- run_imp(a),
                   "1 SNPs are not on the reference genome")
    expect_setequal(res$SNP,
                    c("rs12987662", "rs9320913", "rs11712056", "rs3"))
    expect_equal(res[IMPUTATION_SNP %in% TRUE]$SNP, "rs11712056")
    # b) an indel with the same off-reference ID made the SNV vanish
    b <- rbind(off_ref_dt, data.table::data.table(
        SNP = "rs2000000000", CHR = 3L, BP = 49914397L, A1 = "T", A2 = "TA",
        BETA = 0.1, SE = 0.003, P = 0.5
    ))
    expect_message(res <- run_imp(b),
                   "1 SNPs are not on the reference genome")
    expect_equal(nrow(res), 4)
    expect_equal(res[A2 == "TA"]$SNP, "rs2000000000")
    expect_equal(res[IMPUTATION_SNP %in% TRUE]$SNP, "rs11712056")
})

test_that("#150: NA CHR/BP rows are imputed from the SNP ID", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    skip_if_not_installed("BSgenome.Hsapiens.1000genomes.hs37d5")
    # reporter's data (CHR is character after check_bp_range)
    d <- data.table::data.table(
        SNP = c("rs1000000", "rs9442372"),
        CHR = c("12", NA),
        BP = c(126890980, NA),
        A1 = "G", A2 = "A", FRQ = 0.4, BETA = -0.01, SE = 0.01, P = 1e-8
    )
    run_no_chr_bp <- function(imputation_ind) {
        MungeSumstats:::check_no_chr_bp(
            sumstats_dt = data.table::copy(d), path = "",
            ref_genome = "GRCh37", rsids = NULL,
            imputation_ind = imputation_ind, log_folder_ind = FALSE,
            check_save_out = list(), tabix_index = FALSE, nThread = 1,
            log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
        )$sumstats_dt
    }
    # imputation_ind=TRUE gave an rbindlist column-count error
    res <- suppressMessages(run_no_chr_bp(TRUE))
    expect_equal(nrow(res), 2)
    imputed <- res[SNP == "rs9442372"]
    expect_equal(as.character(imputed$CHR), "1")
    expect_equal(imputed$BP, 1018704)
    expect_true(isTRUE(imputed$IMPUTATION_CHR))
    expect_true(isTRUE(imputed$IMPUTATION_BP))
    expect_true(is.na(res[SNP == "rs1000000"]$IMPUTATION_CHR))
    # imputation_ind=FALSE: same rows, no flag columns (lookup mocked)
    testthat::local_mocked_bindings(
        load_ref_genome_data = function(...) {
            keyed_rsids(SNP = "rs9442372", seqnames = "1", pos = 1018704L)
        },
        .package = "MungeSumstats"
    )
    res <- suppressMessages(run_no_chr_bp(FALSE))
    expect_equal(nrow(res), 2)
    expect_false(any(grepl("^IMPUTATION_", names(res))))
    expect_false(anyNA(res$BP))
})

test_that("#150: check_no_chr_bp accepts a pre-loaded rsids table", {
    d <- data.table::data.table(
        SNP = c("rs1000000", "rs9442372"),
        CHR = c("12", NA),
        BP = c(126890980, NA),
        A1 = "G", A2 = "A", P = 1e-8
    )
    rs <- keyed_rsids(
        SNP = c("rs1000000", "rs9442372"),
        seqnames = c("12", "1"),
        pos = c(126890980L, 1018704L)
    )
    expect_no_error(res <- suppressMessages(MungeSumstats:::check_no_chr_bp(
        sumstats_dt = data.table::copy(d), path = "", ref_genome = "GRCh37",
        rsids = rs, imputation_ind = TRUE, log_folder_ind = FALSE,
        check_save_out = list(), tabix_index = FALSE, nThread = 1,
        log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
    )))
    expect_equal(res$sumstats_dt[SNP == "rs9442372"]$BP, 1018704)
})

test_that("#139: check_no_rs_snp reports how many indels it skips", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    # "foo" (SNV at 1:8490603) and "bar" (indel) both lack an rs ID
    dt <- data.table::data.table(
        SNP = c("rs12987662", "foo", "bar"),
        CHR = c(2L, 1L, 2L),
        BP = c(100821548L, 8490603L, 100821600L),
        A1 = c("A", "T", "A"),
        A2 = c("C", "C", "AT"),
        BETA = c(0.027, 0.019, 0.01),
        SE = 0.003,
        P = c(2.693e-24, 1.794e-08, 0.5)
    )
    msgs <- testthat::capture_messages(
        res <- suppressWarnings(MungeSumstats:::check_no_rs_snp(
            sumstats_dt = data.table::copy(dt), path = "",
            ref_genome = "GRCh37", snp_ids_are_rs_ids = TRUE, indels = TRUE,
            imputation_ind = FALSE, log_folder_ind = FALSE,
            check_save_out = list(), tabix_index = FALSE, nThread = 1,
            log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
        ))$sumstats_dt
    )
    # used to print "Found  Indels" (nrow() of a scalar)
    expect_true(any(grepl("Found 1 Indels", msgs)))
    expect_true(any(grepl("2 SNP IDs are not correctly formatted", msgs)))
    # the indel is dropped rather than given the SNV's rs ID
    expect_setequal(res$SNP, c("rs12987662", "rs301800"))
})

test_that("#177: removed-row logs with NA CHR/BP aren't tabix-sorted", {
    log_dir <- tempfile("logs_")
    cso <- suppressMessages(MungeSumstats:::check_save_path(
        save_path = tempfile(fileext = ".tsv.gz"),
        log_folder = log_dir,
        log_folder_ind = TRUE,
        tabix_index = TRUE
    ))
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", NA), CHR = c("1", "2", NA),
        BP = c(10L, 20L, NA), A1 = "A", A2 = "G", P = c(0.1, 0.2, 0.3)
    )
    expect_no_error(res <- suppressMessages(MungeSumstats:::check_miss_data(
        sumstats_dt = data.table::copy(dt), path = NULL,
        log_folder_ind = TRUE, check_save_out = cso, tabix_index = TRUE,
        nThread = 1, log_files = list(),
        drop_na_cols = c("SNP", "CHR", "BP")
    )))
    expect_equal(nrow(res$sumstats_dt), 2)
    expect_true(file.exists(res$log_files$missing_data))
    logged <- data.table::fread(res$log_files$missing_data)
    expect_equal(nrow(logged), 1)
    expect_true(is.na(logged$BP))
    unlink(log_dir, recursive = TRUE)
})
