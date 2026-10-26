# Regression tests for GitHub issues on inferring the effect allele
# (infer_effect_column checks 1-3). Issues: #168, #186, #188, #193, #208.

# standardise_header() loops over the base mapping plus ~100k effect/frq-allele
# name combinations (~2 s per call). Only rows whose UNCORRECTED name is one of
# the test's headers can ever match, so serving just those rows keeps every
# lookup identical while making each call fast.
local_fast_combns <- function(headers, env = parent.frame()) {
    cmb <- MungeSumstats:::get_eff_frq_allele_combns()
    cmb <- cmb[cmb$UNCORRECTED %in% toupper(headers), ]
    testthat::local_mocked_bindings(
        get_eff_frq_allele_combns = function(...) data.table::copy(cmb),
        .package = "MungeSumstats", .env = env
    )
}

# Check 3 compares A1/A2 with the reference; this serves a 3-SNP reference
# (ref alleles A/C/G, alt alleles G/T/A) instead of loading SNPlocs/BSgenome.
local_fake_ref <- function(env = parent.frame()) {
    ref <- data.table::data.table(
        seqnames = "1", pos = c(100L, 200L, 300L),
        SNP = c("rs1", "rs2", "rs3"), ref_allele = c("A", "C", "G")
    )
    ref[, alt_alleles := list(list("G", "T", "A"))]
    testthat::local_mocked_bindings(
        load_ref_genome_data = function(...) data.table::copy(ref),
        .package = "MungeSumstats", .env = env
    )
}

infer <- function(dt, ...) {
    MungeSumstats::infer_effect_column(
        sumstats_dt = data.table::copy(dt), nThread = 1,
        return_list = FALSE, ...
    )
}

test_that("#168: check 2 infers the effect allele from column names", {
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2"), CHR = 1L, BP = c(100L, 200L),
        A1 = c("G", "T"), A2 = c("A", "C"), BETA1 = c(0.1, -0.2),
        P = c(0.01, 0.5)
    )
    dt_frq <- data.table::copy(dt)
    data.table::setnames(dt_frq, "BETA1", "A1FRQ")
    dt_frq[, BETA := c(0.1, -0.2)]
    dt_b2 <- data.table::copy(dt)
    data.table::setnames(dt_b2, "BETA1", "BETA2")
    local_fast_combns(c(names(dt), names(dt_frq), names(dt_b2)))
    # BETA1 / A1FRQ: effect is on A1, so A1/A2 are swapped (A2 = effect)
    expect_message(res <- infer(dt, on_ref_genome = FALSE),
                   "relate to A1 in the sumstat")
    expect_equal(res$A2, c("G", "T"))
    expect_equal(res$A1, c("A", "C"))
    expect_message(res <- infer(dt_frq, on_ref_genome = FALSE),
                   "relate to A1 in the sumstat")
    expect_equal(res$A2, c("G", "T"))
    # BETA2: effect already on A2, nothing swapped
    expect_message(res <- infer(dt_b2, on_ref_genome = FALSE),
                   "relate to A2 in the sumstat")
    expect_equal(res$A1, c("G", "T"))
    expect_equal(res$A2, c("A", "C"))
    # opt-out leaves the data untouched
    expect_no_message(
        res <- infer(dt, on_ref_genome = FALSE, infer_eff_direction = FALSE)
    )
    expect_equal(res, dt)
    # unambiguous allele headers are never re-inferred
    dt_un <- data.table::copy(dt)
    data.table::setnames(dt_un, c("A1", "A2"),
                         c("NON_EFFECT_ALLELE", "EFFECT_ALLELE"))
    expect_no_message(res <- infer(dt_un, on_ref_genome = FALSE),
                      message = "ambiguous")
    expect_equal(res, dt_un)
    # ambiguous alleles, plain BETA and no reference check: left as is
    dt_amb <- data.table::copy(dt)
    data.table::setnames(dt_amb, "BETA1", "BETA")
    expect_message(res <- infer(dt_amb, on_ref_genome = FALSE),
                   "Can't infer allele columns from sumstats")
    expect_equal(res, dt_amb)
})

test_that("#186: EA/OA allele headers are not treated as ambiguous", {
    input <- data.table::data.table(
        RSID = c("rs1", "rs2"), BETA = c(-0.013, 0.007), SE = c(0.106, 0.10),
        CHISQ = c(0.04, 0.03), P = c(0.77, 0.88), EA = c("A", "T"),
        OA = c("T", "A"), LOG10P = c(0.12, 0.05)
    )
    expect_no_message(res <- infer(input, on_ref_genome = FALSE),
                      message = "ambiguous")
    expect_equal(res, input)
})

test_that("#188/#208: check 3 makes the reference allele the non-effect A1", {
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = 1L, BP = c(100L, 200L, 300L),
        Allele1 = c("A", "C", "G"), Allele2 = c("G", "T", "A"),
        FRQ = c(0.1, 0.2, 0.3), BETA = c(0.1, -0.2, 0.3), SE = 0.01, P = 0.5
    )
    local_fast_combns(names(ss))
    local_fake_ref()
    # Allele1 matches the reference: already A1 = non-effect, no swap
    expect_message(out <- infer(ss, ref_genome = "GRCh37", dbSNP = 144),
                   "relate to A2 in the inputted sumstats")
    expect_equal(names(out), names(ss))
    expect_equal(out$Allele1, ss$Allele1)
    # Allele2 matches the reference: swap so A2 holds the effect allele
    # (#188: the direction of this decision used to be inverted)
    ss_sw <- data.table::copy(ss)
    ss_sw[, c("Allele1", "Allele2") := list(Allele2, Allele1)]
    expect_message(out <- infer(ss_sw, ref_genome = "GRCh37", dbSNP = 144),
                   "relate to A1 in the inputted sumstats")
    expect_equal(out$A2, c("G", "T", "A"))
    expect_equal(out$A1, c("A", "C", "G"))
    # #208: lower-case alleles must still match the reference
    ss_lc <- data.table::copy(ss_sw)
    ss_lc[, c("Allele1", "Allele2") := list(tolower(Allele1),
                                            tolower(Allele2))]
    expect_message(out <- infer(ss_lc, ref_genome = "GRCh37", dbSNP = 144),
                   "relate to A1 in the inputted sumstats")
    expect_equal(out$A2, c("g", "t", "a"))
    expect_equal(out$A1, c("a", "c", "g"))
})

test_that("#193: ALLELE0 headers are ambiguous; A1FREQ marks ALLELE1 effect", {
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = 1L, BP = c(100L, 200L, 300L),
        ALLELE1 = c("A", "C", "G"), ALLELE0 = c("G", "T", "A"),
        A1FREQ = c(0.9, 0.8, 0.7), BETA = c(0.1, 0.2, 0.3), SE = 0.01,
        P = 0.5
    )
    ss_frq <- data.table::copy(ss)
    data.table::setnames(ss_frq, "A1FREQ", "FREQ")
    local_fast_combns(c(names(ss), "FREQ"))
    # A0/A1 naming is ambiguous (pre-fix ALLELE0 was not, so no inference)
    expect_message(res <- infer(ss_frq, ref_genome = "GRCh37",
                                on_ref_genome = FALSE),
                   "Allele columns are ambiguous")
    expect_equal(res, ss_frq)
    # A1FREQ: effect is on ALLELE1, which becomes MSS's A2 (effect)
    expect_message(res <- infer(ss, ref_genome = "GRCh37",
                                on_ref_genome = FALSE),
                   "relate to A1 where A0 in the sumstat")
    expect_equal(names(res), names(ss))
    std <- suppressMessages(
        MungeSumstats::standardise_header(res, return_list = FALSE)
    )
    expect_equal(std$A2, ss$ALLELE1)
    expect_equal(std$A1, ss$ALLELE0)
})

test_that("#193: eff_on_minor_alleles infers from the reference only", {
    ss <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), CHR = 1L, BP = c(100L, 200L, 300L),
        ALLELE1 = c("A", "C", "G"), ALLELE0 = c("G", "T", "A"),
        A1FREQ = c(0.9, 0.8, 0.7), BETA = c(0.1, 0.2, 0.3), SE = 0.01,
        P = 0.5
    )
    expect_error(
        infer(ss, ref_genome = "GRCh37", on_ref_genome = FALSE,
              eff_on_minor_alleles = TRUE),
        "on_ref_genome must equal TRUE to use eff_on_minor_alleles"
    )
    local_fast_combns(names(ss))
    local_fake_ref()
    # ALLELE1 matches the reference, so it is taken as the non-effect A1
    # even though A1FREQ names it (check 2 is skipped)
    expect_message(
        out <- infer(ss, ref_genome = "GRCh37", dbSNP = 144,
                     eff_on_minor_alleles = TRUE),
        "relate to A1 in the inputted sumstats"
    )
    expect_equal(out$A1, c("A", "C", "G"))
    expect_equal(out$A2, c("G", "T", "A"))
})

test_that("#208: mixed-case allele headers are swapped by check 2", {
    x <- data.table::data.table(
        Chr = 1L, Pos_b37 = c(693731L, 723891L),
        RSID = c("rs12238997", "rs2977670"), Allele1 = c("a", "c"),
        Allele2 = c("g", "g"), A1FREQ = c(0.87, 0.86),
        Effect = c(0.0088, -0.0388), StdErr = c(0.0196, 0.0393),
        "P-value" = c(0.6531, 0.3236), n_total_sum = c(334575L, 119680L)
    )
    x_uc <- data.table::copy(x)
    data.table::setnames(x_uc, toupper(names(x_uc)))
    x_a0 <- data.table::data.table(
        SNP = c("rs1", "rs2"), CHR = 1L, BP = c(100L, 200L),
        Allele1 = c("A", "C"), Allele0 = c("G", "T"),
        FRQ_A0 = c(0.9, 0.8), BETA = c(0.1, 0.2), P = 0.5
    )
    local_fast_combns(c(names(x), names(x_a0)))
    # upper-case headers (control) and mixed case (used to error in setnames)
    for (d in list(x_uc, x)) {
        expect_message(out <- infer(d, ref_genome = "GRCh37",
                                    on_ref_genome = FALSE),
                       "relate to A1 in the sumstat")
        expect_equal(out$A2, c("a", "c"))
        expect_equal(out$A1, c("g", "g"))
    }
    # A0 branch with mixed case: the effect allele (Allele0) becomes A2
    expect_message(out <- infer(x_a0, ref_genome = "GRCh37",
                                on_ref_genome = FALSE),
                   "relate to A0 in the sumstat")
    expect_equal(out$A1, c("A", "C"))
    expect_equal(out$A2, c("G", "T"))
})

test_that("#208: raw METAL headers (Freq1) mark Allele1 as the effect", {
    x <- data.table::data.table(
        Chr = 1L, Pos_b37 = c(693731L, 723891L),
        RSID = c("rs12238997", "rs2977670"), Allele1 = c("a", "c"),
        Allele2 = c("g", "g"), Freq1 = c(0.87, 0.86),
        Effect = c(0.0088, -0.0388), StdErr = c(0.0196, 0.0393),
        "P-value" = c(0.6531, 0.3236), n_total_sum = c(334575L, 119680L)
    )
    local_fast_combns(names(x))
    # METAL: Allele1 is the effect allele and Freq1 is its frequency
    expect_message(out <- infer(x, ref_genome = "GRCh37",
                                on_ref_genome = FALSE),
                   "relate to A1 in the sumstat")
    expect_equal(out$A2, c("a", "c"))
    expect_equal(out$A1, c("g", "g"))
})
