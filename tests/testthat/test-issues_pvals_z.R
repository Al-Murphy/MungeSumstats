# Regression tests for GitHub issues on P range, Z-score, multi-trait P and
# BETA imputation checks. Issues: #76, #129, #135, #149, #212, #213.

test_that("#76: out-of-range p-values are set to 1/0 only when converting", {
    dt <- data.table::data.table(SNP = paste0("rs", 1:5),
                                 P = c(5, 0.5, -5, NA, 1))
    run_range <- function(large, neg) {
        MungeSumstats:::check_range_p_val(
            sumstats_dt = data.table::copy(dt), convert_large_p = large,
            convert_neg_p = neg, imputation_ind = TRUE
        )$sumstats_dt
    }
    expect_message(res <- run_range(TRUE, TRUE),
                   "1 p-values are >1 .*These will be converted to 1")
    expect_equal(res$P, c(1, 0.5, 0, NA, 1))
    # opt-out: values untouched and no indicator columns
    expect_message(res <- run_range(FALSE, FALSE),
                   "p-values are <0 .*These will NOT be converted")
    expect_equal(res$P, dt$P)
    expect_false(any(c("convert_large_p_1", "convert_neg_p_0") %in%
                         names(res)))
})

test_that("#212/#213: P=0 gives a finite Z (P set to 1e-300) and is flagged", {
    dt <- data.table::data.table(
        SNP = paste0("rs", 1:3), BETA = c(0.5, -0.2, 0.1),
        SE = c(0.01, 0.02, 0.03), P = c(0, 1e-5, 0.5)
    )
    run_z <- function(d, compute_z) {
        MungeSumstats:::check_zscore(
            sumstats_dt = data.table::copy(d), imputation_ind = TRUE,
            compute_z = compute_z, force_new_z = FALSE,
            standardise_headers = FALSE,
            mapping_file = MungeSumstats:::sumstatsColHeaders
        )$sumstats_dt
    }
    expect_message(res <- run_z(dt, "P"), "Found 1 p-values equal to 0")
    expect_true(all(is.finite(res$Z)))
    expect_equal(res$P[1], 1e-300)
    expect_equal(res$Z[1], sqrt(stats::qchisq(1e-300, 1, lower.tail = FALSE)))
    expect_equal(res$Z[2], -sqrt(stats::qchisq(1e-5, 1, lower.tail = FALSE)))
    expect_true(isTRUE(res$IMPUTATION_p_for_zero[1]))
    expect_true(all(is.na(res$IMPUTATION_p_for_zero[2:3])))
    expect_true(all(res$IMPUTATION_z_score_p))
    # compute_z=TRUE (the reporter's call) takes the same P route
    expect_equal(suppressMessages(run_z(dt, TRUE)), res)
    # the BETA/SE route leaves P=0 alone
    res_b <- run_z(dt, "BETA")
    expect_equal(res_b$Z, dt$BETA / dt$SE)
    expect_equal(res_b$P[1], 0)
    # an NA P must not break the zero check (gives NA Z)
    dt_na <- data.table::copy(dt)[, P := c(0, NA, 0.5)]
    expect_no_error(res <- suppressMessages(run_z(dt_na, "P")))
    expect_true(is.na(res$Z[2]))
    expect_true(is.finite(res$Z[1]))
})

test_that("#149: multi-trait P columns error unless ignore_multi_trait=TRUE", {
    dt <- data.table::data.table(
        SNP = "rs1", CHR = 1L, BP = 100L, A1 = "A", A2 = "G", N = 100L,
        NSTUDY = 3L, FRQ = 0.2, BETA = 0.1, SE = 0.01,
        PVALUE_NEG_LOG10 = 1.30103, P = 0.05, PVALUE_NEG_LOG10_GC = 1.2,
        PVALUE_GC = 0.063
    )
    run_multi <- function(ignore) {
        MungeSumstats:::check_multi_gwas(
            sumstats_dt = data.table::copy(dt), path = NULL,
            analysis_trait = NULL, ignore_multi_trait = ignore,
            mapping_file = MungeSumstats:::sumstatsColHeaders
        )
    }
    expect_error(suppressMessages(run_multi(FALSE)),
                 "ignore_multi_trait=TRUE")
    expect_message(try(run_multi(FALSE), silent = TRUE),
                   "NEG_LOG10, NEG_LOG10_GC, GC")
    # ignoring keeps every column and the plain P value
    expect_no_message(res <- run_multi(TRUE))
    expect_equal(names(res$sumstats_dt), names(dt))
    expect_equal(res$sumstats_dt$P, 0.05)
})

test_that("#129/#135: impute_beta prefers log(OR) > LOG_ODDS > Z*SE, opt-in", {
    d_or <- data.table::data.table(
        SNP = c("rs1", "rs2"), CHR = "1", BP = 1:2, A1 = "A", A2 = "G",
        OR = c(1.5, 0.8), Z = c(2, -1), SE = c(0.1, 0.2), P = c(0.05, 0.3)
    )
    run_signed <- function(d, impute) {
        MungeSumstats:::check_signed_col(
            sumstats_dt = data.table::copy(d), impute_beta = impute,
            log_folder_ind = FALSE, rsids = NULL, imputation_ind = TRUE,
            check_save_out = NULL, tabix_index = FALSE, log_files = list(),
            nThread = 1
        )$sumstats_dt
    }
    # OR wins over LOG_ODDS and Z*SE
    d_both <- data.table::copy(d_or)[, LOG_ODDS := c(9, 9)]
    expect_message(res <- run_signed(d_both, TRUE), "Deriving BETA from OR")
    expect_equal(res$BETA, log(c(1.5, 0.8)))
    expect_true(all(res$IMPUTATION_BETA))
    # LOG_ODDS next
    d_lo <- data.table::copy(d_or)[, OR := NULL][, LOG_ODDS := c(0.4, -0.2)]
    expect_message(res <- run_signed(d_lo, TRUE),
                   "Deriving BETA from LOG ODDS")
    expect_equal(res$BETA, c(0.4, -0.2))
    # Z*SE last
    d_z <- data.table::copy(d_or)[, OR := NULL]
    expect_message(res <- run_signed(d_z, TRUE), "Deriving BETA from Z and SE")
    expect_equal(res$BETA, c(0.2, -0.2))
    # default: no imputation, just a hint
    expect_message(
        res <- run_signed(d_or, FALSE),
        "BETA is not present but can be imputed with OR. Set impute_beta=TRUE"
    )
    expect_false("BETA" %in% names(res))
    # #135: both switches are format_sumstats() arguments
    expect_true(all(c("impute_beta", "es_is_beta") %in%
                        names(formals(MungeSumstats::format_sumstats))))
})
