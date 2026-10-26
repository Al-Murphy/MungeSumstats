# Regression tests for GitHub issues on N, FRQ, INFO and missing-data checks.
# Issues: #3, #5, #26, #29, #43, #58, #100, #151, #180, #181, #182, #189, #197.
cso <- list(log_folder = tempdir(), extension = ".tsv.gz", sep = "\t")

test_that("#26: character N is converted with or without imputation_ind", {
    dt <- data.table::data.table(SNP = c("rs1", "rs2", "rs3"),
                                 N = c("1000", "1500.6", "2000"))
    for (imp in c(FALSE, TRUE)) {
        out <- MungeSumstats:::check_n_int(
            sumstats_dt = data.table::copy(dt), path = NULL,
            convert_n_int = TRUE, imputation_ind = imp
        )$sumstats_dt
        expect_true(is.numeric(out$N))
        expect_equal(out$N, c(1000, 1501, 2000))
    }
    # the imputation flag marks only the rounded value
    expect_equal(out$convert_n_integer, c(NA, TRUE, NA))
})

test_that("#29: N_dropNA=TRUE drops NA N and NA N doesn't disable N_std", {
    run_n_num <- function(dt, drop) {
        MungeSumstats:::check_n_num(
            sumstats_dt = data.table::copy(dt), path = NULL, N_std = 5,
            N_dropNA = drop, log_folder_ind = FALSE, check_save_out = cso,
            tabix_index = FALSE, nThread = 1, log_files = list()
        )$sumstats_dt
    }
    dt <- data.table::data.table(SNP = paste0("rs", 1:3), N = c(100, 100, NA))
    # documented: "N_dropNA Drop rows where N is missing. Default is TRUE."
    expect_equal(nrow(run_n_num(dt, TRUE)), 2)
    expect_equal(nrow(run_n_num(dt, FALSE)), 3)
    # one NA N must not switch off the outlier (N_std) filter
    dt2 <- data.table::data.table(SNP = paste0("rs", 1:41),
                                  N = c(rep(100, 39), 1e6, NA))
    out <- run_n_num(dt2, TRUE)
    expect_false("rs40" %in% out$SNP)
    expect_equal(nrow(out), 39)
})

test_that("#43: compute_nsize formulas, Neff_<method> names, IMPUTATION_n", {
    # powers of two keep the arithmetic exact
    dtA <- data.table::data.table(SNP = c("rs1", "rs2"),
                                  N_CAS = c(1024, 1024), N_CON = c(1024, 4096))
    dtB <- data.table::data.table(SNP = c("rs1", "rs2"),
                                  N_CAS = c(1024, 1024), N_CON = c(1024, 3072))
    nsize <- function(dt, m, ...) {
        suppressMessages(MungeSumstats::compute_nsize(
            sumstats_dt = dt, compute_n = m, return_list = FALSE, ...
        ))
    }
    # giant: 2/(1/N_CAS + 1/N_CON); metal: 4/(1/N_CAS + 1/N_CON)
    expect_equal(nsize(dtA, "giant")$Neff, c(1024L, 1638L))
    expect_equal(nsize(dtA, "metal")$Neff, c(2048L, 3276L))
    expect_equal(nsize(dtA, "sum")$N, c(2048L, 5120L))
    # ldsc: (N_CAS+N_CON) * P / mean(P at max(N_CAS+N_CON))
    expect_equal(nsize(dtB, "ldsc")$Neff, c(4096L, 4096L))
    multi <- nsize(dtA, c("ldsc", "giant", "metal"))
    expect_true(all(c("Neff_ldsc", "Neff_giant", "Neff_metal") %in%
                        names(multi)))
    expect_false("Neff" %in% names(multi))
    expect_equal(multi$Neff_giant, c(1024L, 1638L))
    # numeric compute_n sets N and, if asked, the imputation flag
    out <- nsize(data.table::data.table(SNP = "rs1", P = 0.1), 1000,
                 imputation_ind = TRUE)
    expect_equal(out$N, 1000)
    expect_true("IMPUTATION_n" %in% names(out))
    # an existing N is kept and not flagged as imputed
    out <- nsize(data.table::data.table(SNP = "rs1", P = 0.1, N = 5), 1000,
                 imputation_ind = TRUE)
    expect_equal(out$N, 5)
    expect_false("IMPUTATION_n" %in% names(out))
})

test_that("#3/#151: check_ldsc_format forces the parameters LDSC needs", {
    dt <- data.table::data.table(SNP = "rs1", A1 = "A", A2 = "G", BETA = 0.1,
                                 SE = 0.01, P = 0.01, N = 1000)
    ldsc <- function(d, save_format = "LDSC", compute_n = 0L) {
        suppressMessages(MungeSumstats:::check_ldsc_format(
            sumstats_dt = d, save_format = save_format,
            convert_n_int = FALSE, allele_flip_check = FALSE,
            compute_z = FALSE, compute_n = compute_n
        ))
    }
    for (sf in c("LDSC", "ldsc")) {
        out <- ldsc(dt, save_format = sf)
        expect_true(out$convert_n_int)
        expect_true(out$allele_flip_check)
        expect_equal(out$compute_z, "BETA")
    }
    # each forced parameter is reported
    expect_message(
        MungeSumstats:::check_ldsc_format(
            sumstats_dt = dt, save_format = "LDSC", convert_n_int = TRUE,
            allele_flip_check = FALSE, compute_z = "BETA", compute_n = 0L
        ),
        "Setting `allele_flip_check=TRUE`", fixed = TRUE
    )
    # Z from P when BETA/SE are missing
    expect_equal(ldsc(dt[, !c("BETA", "SE")])$compute_z, "P")
    # LDSC needs N: missing N errors unless compute_n will add it
    expect_error(ldsc(dt[, !"N"]), "LDSC requires an N column")
    expect_equal(ldsc(dt[, !"N"], compute_n = 1000)$compute_z, "BETA")
    # other formats are left alone
    out <- ldsc(dt, save_format = NULL)
    expect_false(out$convert_n_int)
    expect_false(out$allele_flip_check)
    expect_false(out$compute_z)
})

test_that("#3: LDSC output accepts several compute_n methods", {
    dt <- data.table::data.table(SNP = "rs1", A1 = "A", A2 = "G", BETA = 0.1,
                                 SE = 0.01, P = 0.01, N_CAS = 10, N_CON = 20)
    expect_no_error(suppressMessages(MungeSumstats:::check_ldsc_format(
        sumstats_dt = dt, save_format = "LDSC", convert_n_int = FALSE,
        allele_flip_check = FALSE, compute_z = FALSE,
        compute_n = c("ldsc", "giant")
    )))
})

test_that("#58/#189: check_frq_maf renames or warns, never converts FRQ", {
    dt <- data.table::data.table(SNP = c("rs1", "rs2"), FRQ = c(0.7, 0.2))
    # frq_is_maf=FALSE: values > 0.5 mean the column is a major-allele FRQ
    out <- suppressMessages(MungeSumstats:::check_frq_maf(
        sumstats_dt = data.table::copy(dt), frq_is_maf = FALSE
    ))
    expect_true("MAJOR_ALLELE_FRQ" %in% names(out))
    expect_false("FRQ" %in% names(out))
    expect_equal(out$MAJOR_ALLELE_FRQ, dt$FRQ)
    # frq_is_maf=TRUE: warn only, no 1-FRQ conversion
    expect_message(
        out <- MungeSumstats:::check_frq_maf(
            sumstats_dt = data.table::copy(dt), frq_is_maf = TRUE
        ),
        "1 SNPs \\(50%\\) have FRQ values > 0.5"
    )
    expect_equal(out$FRQ, dt$FRQ)
    # nothing > 0.5: FRQ kept even with frq_is_maf=FALSE
    out <- MungeSumstats:::check_frq_maf(
        sumstats_dt = data.table::data.table(SNP = "rs1", FRQ = 0.2),
        frq_is_maf = FALSE
    )
    expect_equal(names(out), c("SNP", "FRQ"))
})

test_that("#180/#181/#182/#197: default drop_na_cols ignores extra columns", {
    dt <- data.table::data.table(
        SNP = paste0("rs", 1:4), CHR = c("1", "1", "2", "2"),
        BP = c(100L, 200L, 300L, 400L), A1 = c("A", "C", "G", "T"),
        A2 = c("G", "T", "A", "C"), BETA = c(0.1, -0.2, 0.3, NA),
        P = c(0.01, 0.02, 0.03, 0.04), EXTRA = NA,
        ROW_NUMBER = c(1L, NA, 3L, 4L), ANNOTATION = c("a", NA, "c", "d")
    )
    miss <- function(cols) {
        MungeSumstats:::check_miss_data(
            sumstats_dt = data.table::copy(dt), path = NULL,
            log_folder_ind = FALSE, check_save_out = NULL,
            tabix_index = FALSE, nThread = 1, log_files = list(),
            drop_na_cols = cols
        )$sumstats_dt
    }
    # the format_sumstats() default (check_miss_data() has none)
    default_cols <- eval(formals(MungeSumstats::format_sumstats)$drop_na_cols)
    expect_message(
        out <- miss(default_cols),
        "WARNING: 1 rows in sumstats file are missing data and will be removed"
    )
    # only the NA-BETA row goes; extra columns and rs2's NA extras are kept
    expect_equal(out$SNP, c("rs1", "rs2", "rs3"))
    expect_true(all(c("EXTRA", "ROW_NUMBER", "ANNOTATION") %in% names(out)))
    # NULL checks every column, so the all-NA EXTRA removes every row
    expect_error(suppressMessages(miss(NULL)),
                 "All SNPs have been filtered out")
    # a custom set is used as given
    out <- suppressMessages(miss("ROW_NUMBER"))
    expect_equal(out$SNP, c("rs1", "rs3", "rs4"))
})

test_that("#5/#100: an all-NA or absent INFO column drops no rows", {
    fe <- data.table::data.table(SNP = paste0("rs", 1:4), CHR = "1",
                                 BP = 1:4, P = 0.01)
    run_info <- function(dt) {
        MungeSumstats:::check_info_score(
            sumstats_dt = dt, INFO_filter = 0.9, log_folder_ind = TRUE,
            check_save_out = cso, tabix_index = FALSE, nThread = 1,
            log_files = list()
        )
    }
    na_info <- data.table::copy(fe)[, INFO := NA_real_]
    out <- suppressMessages(run_info(na_info))
    expect_equal(nrow(out$sumstats_dt), nrow(fe))
    expect_false("info_filter" %in% names(out$log_files))
    # #100: no INFO column (e.g. a VCF whose INFO only holds AF) skips the step
    expect_message(out <- run_info(data.table::copy(fe)),
                   "INFO column not available")
    expect_equal(nrow(out$sumstats_dt), nrow(fe))
    expect_length(out$log_files, 0)
})
