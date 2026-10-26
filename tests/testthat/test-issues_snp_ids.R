# Regression tests for GitHub issues on SNP ID parsing.
# Issues: #120, #123 (non-rs SNP IDs), #138 (combined CHR:BP:A1:A2 IDs),
# #48, #57 (multi-rs and CHR:BP columns), #159/#162/#176 (imputation_ind).
run_check_no_rs_snp <- function(dt) {
    # non-rs IDs raise "NAs introduced by coercion" warnings by design
    suppressWarnings(MungeSumstats:::check_no_rs_snp(
        sumstats_dt = data.table::copy(dt),
        path = "",
        ref_genome = "GRCh37",
        snp_ids_are_rs_ids = TRUE,
        indels = TRUE,
        imputation_ind = FALSE,
        log_folder_ind = FALSE,
        check_save_out = NULL,
        tabix_index = FALSE,
        nThread = 1,
        log_files = list(),
        dbSNP = 144,
        dbSNP_tarball = NULL
    ))
}

test_that("#120: mix of rs IDs and CHR:BP:A1:A2 IDs doesn't crash", {
    dt <- data.table::data.table(
        SNP = c("chr1:701203:G:T", "rs185127847", "chr1:722408:G:C"),
        A1 = c("T", "A", "C"),
        A2 = c("G", "T", "G"),
        BETA = c(1.2, 0.001, 0.01),
        P = c(0.8, 0.65, 0.45)
    )
    expect_message(
        res <- run_check_no_rs_snp(dt),
        "2 SNP IDs are not correctly formatted and will be removed"
    )
    expect_equal(res$sumstats_dt$SNP, "rs185127847")
})

test_that("#123: all CHR:BP SNP IDs are moved to ID, not looked up", {
    dt <- data.table::data.table(
        SNP = c("4:32435284", "16:17044975", "1:53778300"),
        CHR = c(4L, 16L, 1L),
        BP = c(32435284L, 17044975L, 53778300L),
        A1 = c("A", "T", "A"),
        A2 = c("G", "C", "G"),
        P = c(1.673e-09, 6.265e-09, 1.527e-08)
    )
    expect_no_message(
        res <- run_check_no_rs_snp(dt),
        message = "Checking SNP RSIDs"
    )
    expect_true("ID" %in% names(res$sumstats_dt))
    expect_false("SNP" %in% names(res$sumstats_dt))
    expect_equal(res$sumstats_dt$ID, dt$SNP)
})

test_that("#138: combined CHR:BP:A?:A? column is split in the right order", {
    run_four_step <- function(dt) {
        MungeSumstats:::check_four_step_col(
            sumstats_dt = data.table::copy(dt),
            path = NULL
        )$sumstats_dt
    }
    # a) no allele columns: default order is CHR:BP:A1:A2
    a <- data.table::data.table(
        ID = c("1:15791:C:T", "1:69487:G:A"),
        BETA = c(-3.56, -0.011)
    )
    expect_message(
        res <- run_four_step(a),
        "separated into the columns CHR, BP, A1, A2"
    )
    expect_equal(res$A1, c("C", "G"))
    expect_equal(res$A2, c("T", "A"))
    expect_equal(as.integer(res$BP), c(15791L, 69487L))
    expect_false("ID" %in% names(res))
    # b) separate A2 column matches the 4th field -> CHR:BP:A1:A2
    b <- data.table::data.table(
        ID = c("1:15791:C:T", "1:69487:G:A"),
        A2 = c("T", "A"),
        BETA = c(-3.56, -0.011)
    )
    res <- run_four_step(b)
    expect_equal(res$A1, c("C", "G"))
    expect_equal(res$A2, c("T", "A"))
    # c) separate A2 column matches the 3rd field -> CHR:BP:A2:A1
    c_dt <- data.table::data.table(
        ID = c("1:15791:T:C", "1:69487:A:G"),
        A2 = c("T", "A")
    )
    res <- run_four_step(c_dt)
    expect_equal(res$A1, c("C", "G"))
    expect_equal(res$A2, c("T", "A"))
    # d) order declared in the column name
    d <- data.table::data.table(`CHR:BP:A2:A1` = c("1:15791:T:C",
                                                   "1:69487:A:G"))
    res <- run_four_step(d)
    expect_equal(res$A1, c("C", "G"))
    expect_equal(res$A2, c("T", "A"))
})

test_that("#48: multi-rs IDs split on ',', '; ' and '_' separators", {
    dt <- data.table::data.table(
        SNP = c("rs1,rs2", "rs3; rs4", "rs5_rs6", "rs7"),
        P = 0.1
    )
    run_multi <- function(remove) {
        MungeSumstats:::check_multi_rs_snp(
            sumstats_dt = data.table::copy(dt), path = NULL,
            remove_multi_rs_snp = remove, imputation_ind = TRUE,
            log_folder_ind = FALSE, check_save_out = list(),
            tabix_index = FALSE, nThread = 1, log_files = list()
        )$sumstats_dt
    }
    # FALSE: keep the first rs ID of each row and flag it
    expect_message(
        res <- run_multi(FALSE),
        "3 SNPs found with multiple RSIDs on one row, the first will be taken"
    )
    expect_equal(res$SNP, c("rs1", "rs3", "rs5", "rs7"))
    expect_equal(res$convert_multi_rs_SNP, c(TRUE, TRUE, TRUE, NA))
    # TRUE: drop every multi-rs row whatever the separator
    res <- suppressMessages(run_multi(TRUE))
    expect_equal(res$SNP, "rs7")
    expect_false("convert_multi_rs_SNP" %in% names(res))
})

test_that("#48: remove_multi_rs_snp=TRUE message counts all separators", {
    dt <- data.table::data.table(
        SNP = c("rs1,rs2", "rs3; rs4", "rs5_rs6", "rs7"),
        P = 0.1
    )
    expect_message(
        MungeSumstats:::check_multi_rs_snp(
            sumstats_dt = data.table::copy(dt), path = NULL,
            remove_multi_rs_snp = TRUE, imputation_ind = FALSE,
            log_folder_ind = FALSE, check_save_out = list(),
            tabix_index = FALSE, nThread = 1, log_files = list()
        ),
        "^3 SNPs found with multiple RSIDs on one row, these will be removed"
    )
})

test_that("#48: CHR:BP column with a 3-field CHR:BP:ID row keeps CHR, BP", {
    dt <- data.table::data.table(
        MarkerName = c("11:107819621", "5:85928892", "16:23609681:ID"),
        P = c(0.1, 0.2, 0.3)
    )
    expect_message(
        res <- MungeSumstats:::check_two_step_col(
            sumstats_dt = data.table::copy(dt),
            path = NULL
        )$sumstats_dt,
        "MarkerName has been separated into the columns CHR, BP"
    )
    expect_equal(names(res), c("P", "CHR", "BP"))
    expect_equal(res$CHR, c(11L, 5L, 16L))
    expect_equal(res$BP, c(107819621L, 85928892L, 23609681L))
})

test_that("#57: VARIANT CHR:BP:A1:A2 column gives integer BP for joins", {
    dt <- data.table::data.table(
        VARIANT = c("1:69487:G:A", "10:69569:T:C"),
        BETA = c(0.1, -0.2),
        P = c(0.5, 0.6)
    )
    res <- suppressMessages(MungeSumstats:::check_four_step_col(
        sumstats_dt = data.table::copy(dt),
        path = NULL
    ))$sumstats_dt
    expect_false("VARIANT" %in% names(res))
    expect_true(is.integer(res$CHR))
    expect_true(is.integer(res$BP))
    expect_equal(res$BP, c(69487L, 69569L))
    expect_equal(res$A1, c("G", "T"))
    expect_equal(res$A2, c("A", "C"))
    # the reported "Incompatible join types: x.BP (character) and
    # i.BP (integer)" came from joining this against integer positions
    ref <- data.table::data.table(CHR = 10L, BP = 69569L, hit = TRUE)
    expect_no_error(joined <- ref[res, on = c("CHR", "BP")])
    expect_equal(joined$hit, c(NA, TRUE))
})

test_that("#176: NA SNP IDs with imputation_ind=TRUE get unique join keys", {
    # stand-in for the reference lookup: every row gets a new rs ID
    testthat::local_mocked_bindings(
        check_no_snp = function(sumstats_dt, log_files, ...) {
            sumstats_dt[, SNP := paste0("rs", 9000 + seq_len(.N))]
            sumstats_dt[, IMPUTATION_SNP := TRUE]
            list(sumstats_dt = sumstats_dt, log_files = log_files)
        },
        .package = "MungeSumstats"
    )
    # 5 NA keys joined to 5 NA keys = 25 rows > 6 + 5, which raised
    # data.table's cartesian-join error before the fix
    dt <- data.table::data.table(
        SNP = c("rs1", NA, "", NA, NA, NA),
        CHR = "1",
        BP = seq(100L, 600L, 100L),
        A1 = "A",
        A2 = "G",
        P = 0.5
    )
    expect_message(
        res <- MungeSumstats:::check_no_rs_snp(
            sumstats_dt = data.table::copy(dt), path = NULL,
            ref_genome = "GRCh37", snp_ids_are_rs_ids = TRUE, indels = TRUE,
            imputation_ind = TRUE, log_folder_ind = FALSE,
            check_save_out = NULL, tabix_index = FALSE, nThread = 1,
            log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
        )$sumstats_dt,
        "5 SNP IDs are not correctly formatted"
    )
    expect_equal(nrow(res), 6)
    expect_true(all(grepl("^rs", res$SNP)))
    expect_false("SNP_old_temp" %in% names(res))
    expect_equal(sum(res$IMPUTATION_SNP, na.rm = TRUE), 5)
    expect_true(is.na(res[SNP == "rs1"]$IMPUTATION_SNP))
})

test_that("#159/#162: imputation_ind=TRUE with chr:bp and malformed IDs", {
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    run_imp <- function(dt) {
        # non-rs IDs raise "NAs introduced by coercion" warnings by design
        suppressWarnings(MungeSumstats:::check_no_rs_snp(
            sumstats_dt = data.table::copy(dt), path = "",
            ref_genome = "GRCh37", snp_ids_are_rs_ids = TRUE, indels = TRUE,
            imputation_ind = TRUE, log_folder_ind = FALSE,
            check_save_out = list(), tabix_index = FALSE, nThread = 1,
            log_files = list(), dbSNP = 144, dbSNP_tarball = NULL
        ))$sumstats_dt
    }
    dt <- data.table::data.table(
        SNP = c("rs12987662", "rs9320913", "3:49914397", "foo"),
        CHR = c(2L, 6L, 3L, 1L),
        BP = c(100821548L, 98584733L, 49914397L, 8490603L),
        A1 = c("A", "A", "T", "T"),
        A2 = c("C", "C", "C", "C"),
        BETA = c(0.027, 0.024, 0.024, 0.019),
        SE = 0.003,
        P = c(2.693e-24, 2.457e-19, 3.304e-19, 1.794e-08)
    )
    check_res <- function(res) {
        expect_equal(nrow(res), 4)
        expect_setequal(
            res$SNP,
            c("rs12987662", "rs9320913", "rs11712056", "rs301800")
        )
        expect_false("SNP_old_temp" %in% names(res))
        expect_setequal(res[IMPUTATION_SNP %in% TRUE]$SNP,
                        c("rs11712056", "rs301800"))
        expect_true(all(is.na(
            res[SNP %in% c("rs12987662", "rs9320913")]$IMPUTATION_SNP
        )))
    }
    # #159: rs + chr:bp + malformed ID (crashed on missing IMPUTATION_SNP)
    msgs <- testthat::capture_messages(res <- run_imp(dt))
    expect_true(any(grepl("1 SNP IDs are not correctly formatted", msgs)))
    expect_true(any(grepl("1 SNP IDs appear to be made up of chr:bp", msgs)))
    check_res(res)
    # #162: rs + chr:bp IDs only (rbindlist column-count mismatch)
    dt[4, SNP := "1:8490603"]
    expect_message(res <- run_imp(dt),
                   "2 SNP IDs appear to be made up of chr:bp")
    check_res(res)
})
