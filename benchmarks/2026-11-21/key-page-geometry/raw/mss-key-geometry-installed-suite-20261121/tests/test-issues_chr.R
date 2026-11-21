# Regression tests for GitHub issues on chromosome naming/standardisation.
# Issues: #103, #116, #130, #31, #37, #157, #158 (CHR standardisation),
# #8, #116/#117 (sorting), #144, #150, #161 (BP range and column types).
test_that("#103/#116/#130: check_chr standardises X, Y and MT names", {
    dt <- data.table::data.table(
        SNP = paste0("rs", 1:9),
        CHR = c("1", "x", "chrX", "23", "Y", "24", "MT", "chrM", "chrMT"),
        BP = seq(100L, 900L, 100L)
    )
    # rmv_chr=NULL keeps the sex and mitochondrial chromosomes
    expect_no_message(
        res <- MungeSumstats:::check_chr(
            sumstats_dt = data.table::copy(dt),
            rmv_chr = NULL
        ),
        message = "nonstandard CHR"
    )
    expect_equal(nrow(res$sumstats_dt), 9)
    expect_equal(
        res$sumstats_dt$CHR,
        c("1", "X", "X", "X", "Y", "Y", "MT", "MT", "MT")
    )
    # default style removal drops all of X, Y and MT
    expect_message(
        res2 <- MungeSumstats:::check_chr(
            sumstats_dt = data.table::copy(dt),
            rmv_chr = c("X", "Y", "MT")
        ),
        "8 SNPs are on chromosomes X, Y, MT and will be removed"
    )
    expect_equal(res2$sumstats_dt$CHR, "1")
})

test_that("#116/#117: sort_coords keeps MT and sorts CHR naturally", {
    dt <- data.table::data.table(
        SNP = c("rs2853502", "rs1", "rs2", "rs3"),
        CHR = c("MT", "1", "X", "10"),
        BP = c(13276, 300, 200, 100)
    )
    out <- MungeSumstats:::sort_coords(data.table::copy(dt))
    expect_false(anyNA(out$CHR))
    expect_equal(as.character(out$CHR), c("1", "10", "X", "MT"))
    expect_true(is.integer(out$BP))
})

test_that("#8: sort_coords orders CHR naturally, sex chromosomes last", {
    dt <- data.table::data.table(
        SNP = paste0("rs", 1:6),
        CHR = c("x", "10", "2", "1", "Y", "2"),
        BP = c(5, 4, 30, 2, 1, 3)
    )
    for (method in c("data.table", "GenomicRanges")) {
        out <- suppressMessages(MungeSumstats:::sort_coords(
            sumstats_dt = data.table::copy(dt),
            sort_method = method
        ))
        # lexical order would give 1, 10, 2, 2, X, Y
        expect_equal(out$SNP, c("rs4", "rs6", "rs3", "rs2", "rs1", "rs5"))
        expect_equal(out$CHR, c("1", "2", "2", "10", "X", "Y"))
        expect_true(is.character(out$CHR))
        expect_true(is.integer(out$BP))
    }
})

test_that("#31/#37/#157/#158: check_chr prefixes, contigs and lowercase", {
    dt <- data.table::data.table(
        SNP = paste0("rs", 1:12),
        CHR = c("chr1", "CHR2", "Chr3", "ch4", "22", "chr6_random",
                "GL000192.1", "chr22_KI270731v1_random", "y", "mt", "M",
                "chrx"),
        BP = 1:12
    )
    expect_message(
        res <- MungeSumstats:::check_chr(
            sumstats_dt = data.table::copy(dt),
            rmv_chr = NULL
        ),
        "Removing 3 SNPs with nonstandard CHR entries."
    )
    expect_equal(res$sumstats_dt$SNP, paste0("rs", c(1:5, 9:12)))
    expect_equal(
        res$sumstats_dt$CHR,
        c("1", "2", "3", "4", "22", "Y", "MT", "MT", "X")
    )
    # rmv_chr = c("X", "Y") keeps MT (#158: MT used to become "MTT")
    res2 <- suppressMessages(MungeSumstats:::check_chr(
        sumstats_dt = data.table::copy(dt),
        rmv_chr = c("X", "Y")
    ))
    expect_equal(res2$sumstats_dt$CHR, c("1", "2", "3", "4", "22", "MT", "MT"))
})

run_bp_range <- function(dt, ref_genome = "GRCh37") {
    MungeSumstats:::check_bp_range(
        sumstats_dt = data.table::copy(dt), path = "",
        ref_genome = ref_genome, log_folder_ind = FALSE,
        imputation_ind = FALSE, check_save_out = list(),
        tabix_index = FALSE, nThread = 1, log_files = list()
    )$sumstats_dt
}

test_that("#144: check_bp_range drops BP outside 1..chromosome length", {
    # chr1 length: GRCh37 249,250,621; GRCh38 248,956,422
    dt <- data.table::data.table(
        SNP = paste0("rs", 1:8),
        CHR = c("1", "1", "1", "chr1", "Chr1", NA, "1", "1"),
        BP = c(1000, 2.24e14, 0, 2.24e14, -5, 5, 249000000, NA),
        A1 = "A", A2 = "G", BETA = 0.1, P = 0.5
    )
    expect_message(
        res <- run_bp_range(dt),
        paste0("5 SNPs have been removed as their BP column is not in ",
               "the range of 1 to the length of the chromosome")
    )
    # rs4/rs5: "chr"-prefixed CHR are range-checked too (crash in #144);
    # rs8: NA BP on a known CHR is removed
    # rs6: unknown CHR is kept for check_no_chr_bp (#150)
    expect_equal(res$SNP, c("rs1", "rs6", "rs7"))
    expect_false("chr_max" %in% names(res))
    res38 <- suppressMessages(run_bp_range(dt, ref_genome = "GRCh38"))
    expect_equal(res38$SNP, c("rs1", "rs6"))
})

test_that("#144: check_bp_range also checks 23/24/x/y/M/mt CHR names", {
    chrs <- c("23", "24", "x", "y", "M", "mt", "CHR23", "23")
    dt <- data.table::data.table(
        SNP = paste0("rs", seq_along(chrs)),
        CHR = chrs,
        BP = c(rep(2.24e14, 7), 1000)
    )
    expect_message(
        res <- run_bp_range(dt),
        "^7 SNPs have been removed"
    )
    expect_equal(res$SNP, "rs8")
})

test_that("#144/#150: NA CHR and NA BP rows survive check_bp_range", {
    # used to be dropped silently (not counted or logged) whenever another
    # row had a bad BP, so check_no_chr_bp could never impute it
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs9442372"),
        CHR = c("1", "1", NA),
        BP = c(1000, 2.24e14, NA)
    )
    expect_message(
        res <- run_bp_range(dt),
        "^1 SNPs have been removed"
    )
    expect_equal(res$SNP, c("rs1", "rs9442372"))
})

test_that("#161: character BP and P columns are compared numerically", {
    # all-character columns, as when a header line is repeated in the file
    dt <- data.table::data.table(
        CHR = c("1", "1", "1", "1", "X"),
        BP = c("1", "1000", "2000", "8490603", "X"),
        A1 = c("C", "G", "G", "T", "X"),
        A2 = c("G", "C", "C", "C", "X"),
        P = c("0.05", "0.05", "0.05", "0.05", "X")
    )
    # "2000" > "249250621" as strings; only the "X" row should go
    expect_message(
        res <- suppressWarnings(run_bp_range(dt)),
        "Coercing BP column to numeric."
    )
    expect_true(is.numeric(res$BP))
    expect_equal(res$BP, c(1, 1000, 2000, 8490603))
    # "0.05" < 5e-8 is TRUE as strings, so 4 rows were counted significant
    expect_message(
        suppressWarnings(MungeSumstats:::report_summary(
            sumstats_dt = data.table::copy(dt)
        )),
        "- 0 genome-wide significant variants"
    )
})
