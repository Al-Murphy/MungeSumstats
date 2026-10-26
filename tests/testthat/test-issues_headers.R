# Regression tests for GitHub issues on column-header mapping (standardise_header / sumstatsColHeaders).
# Issues: #12, #45, #48, #58, #91, #94-#97, #99, #105, #138, #142, #148, #168, #174, #175, #186, #192, #195

# standardise_header() appends ~100k allele-specific effect/FRQ name combinations
# (get_eff_frq_allele_combns()) to the map and loops over every row, ~2 s per call.
# A row whose UNCORRECTED never equals a column name (an input name, or a standard
# name a column was renamed to) can never fire. So trimming the combinations to those
# names gives an identical result in milliseconds. The #142 block checks this against
# an untrimmed call. If the refactor removes the internal, std() runs the real thing.
hdr_combns <- if (exists("get_eff_frq_allele_combns",
                         envir = asNamespace("MungeSumstats"),
                         inherits = FALSE)) {
    MungeSumstats:::get_eff_frq_allele_combns()
}
std <- function(dt, ...) {
    run <- function() {
        suppressMessages(MungeSumstats::standardise_header(
            sumstats_dt = data.table::copy(dt), return_list = FALSE, ...
        ))
    }
    if (is.null(hdr_combns)) return(run())
    keep <- c(toupper(names(dt)), as.character(unlist(list(...)$mapping_file)),
              MungeSumstats:::sumstatsColHeaders$Corrected,
              hdr_combns$CORRECTED)
    testthat::with_mocked_bindings(
        run(),
        get_eff_frq_allele_combns = function(...) {
            hdr_combns[hdr_combns$UNCORRECTED %in% keep, ]
        },
        .package = "MungeSumstats"
    )
}

test_that("#12/#148: sumstatsColHeaders ships once, unique and upper-case", {
    m <- MungeSumstats:::sumstatsColHeaders
    expect_true(is.data.frame(m))
    expect_equal(names(m), c("Uncorrected", "Corrected"))
    expect_true(nrow(m) > 400)
    # the internal copy (sysdata.rda) and data("sumstatsColHeaders") must not drift
    e <- new.env()
    utils::data("sumstatsColHeaders", package = "MungeSumstats", envir = e)
    expect_identical(e$sumstatsColHeaders, m)
    # each raw header maps to exactly one standard name; standardise_header()
    # upper-cases the input names, so a lower-case entry could never match
    expect_equal(anyDuplicated(m$Uncorrected), 0)
    expect_equal(m$Uncorrected, toupper(m$Uncorrected))
    expect_equal(m$Uncorrected, trimws(m$Uncorrected))
    expect_false(anyNA(m) || any(m$Uncorrected == "") || any(m$Corrected == ""))
    # the allele-specific combinations never contradict themselves or the base map
    cmb <- MungeSumstats:::get_eff_frq_allele_combns()
    pairs <- unique(data.frame(u = cmb$UNCORRECTED, c = cmb$CORRECTED))
    expect_equal(anyDuplicated(pairs$u), 0)
    shared <- intersect(cmb$UNCORRECTED, m$Uncorrected)
    expect_equal(cmb$CORRECTED[match(shared, cmb$UNCORRECTED)],
                 m$Corrected[match(shared, m$Uncorrected)])
})

test_that("#58/#91/#94/#96/#97/#99/#105/#138/#142/#148/#168/#186: single headers map to the agreed name", {
    # Every mapping decision made in these threads. An expected value equal to
    # toupper(raw) means the header is deliberately left unmapped.
    dec <- do.call(rbind, list(
        data.frame(issue = "#91/#94/#96", expected = "CHR",
                   raw = c("chromosomes", "chro", "seq", "seqnames", "seqs",
                           "seqname", "chroms", "#CHROM")),
        # START is a position; END (second coordinate) is not mapped yet
        data.frame(issue = "#96", raw = c("start", "end"),
                   expected = c("BP", "END")),
        # NGT (number of genotyped cohorts) is not N
        data.frame(issue = "#97", raw = "NGT", expected = "NGT"),
        data.frame(issue = "#58", raw = "AF", expected = "FRQ"),
        data.frame(issue = "#99",
                   raw = c("PosGRCh37", "PosGRCh38", "testedAllele",
                           "otherAllele"),
                   expected = c("BP", "BP", "A2", "A1")),
        # FRQ_U (controls) is the EAF; FRQ_A (cases) stays as is
        data.frame(issue = "#105", raw = c("FRQ_U", "FRQ_A"),
                   expected = c("FRQ", "FRQ_A")),
        data.frame(issue = "#138",
                   raw = c("variant", "minor_allele", "minor_AF"),
                   expected = c("SNP", "A2", "FRQ")),
        data.frame(issue = "#142", expected = "BP",
                   raw = c("GENPOS", "GenomePos", "POSGEN", "GEN_POS",
                           "GEN.POS", "GEN-POS")),
        data.frame(issue = "#148", expected = "A1",
                   raw = c("MAJORALLELE", "major_allele", "MAJOR-ALLELE",
                           "Major Allele")),
        data.frame(issue = "#168",
                   raw = c("TOTALSAMPLESIZE", "VARIANT_ID", "RSIDS",
                           "P_BOLT_LMM", "NCASES", "N_EFFECTIVE",
                           "N_INFORMATIVE", "TOTAL_N", "HET_P", "ALL_AF",
                           "DIRECT", "ALT_EFFSIZE", "INFORMATIVE_ALT_AC"),
                   expected = c("N", "SNP", "SNP", "P", "N_CAS", "N", "N",
                                "N", "HETPVAL", "FRQ", "DIRECTION", "BETA",
                                "AC")),
        # removed in 1.9.18: ambiguous (NMISS, WEIGHT) or GRCh36 coordinates
        data.frame(issue = "#168",
                   raw = c("NMISS", "WEIGHT", "CHR36", "POS36", "BASE_GRCH36"),
                   expected = c("NMISS", "WEIGHT", "CHR36", "POS36",
                                "BASE_GRCH36")),
        # allele-specific names still map, via get_eff_frq_allele_combns()
        data.frame(issue = "#168",
                   raw = c("BETA1", "A1FRQ", "OR_A1", "FREQ.A1.1000G.EUR"),
                   expected = c("BETA", "FRQ", "OR", "FRQ")),
        data.frame(issue = "#186", raw = c("EA", "OA"),
                   expected = c("A2", "A1"))
    ))
    got <- vapply(dec$raw, function(h) {
        dt <- data.table::data.table(x = 1)
        data.table::setnames(dt, h)
        names(std(dt))
    }, character(1))
    names(got) <- paste(dec$issue, dec$raw)
    expect_equal(got, stats::setNames(dec$expected, names(got)))
})

test_that("#95/#105/#138: the first synonym in map order wins and names stay unique", {
    # standardise_header() skips a synonym once its target name exists. This
    # relies on `column_headers` aliasing the data.table's names vector, which
    # setnames() updates in place. Snapshotting it with copy() would rename both
    # AF and MINOR_AF to FRQ (duplicate columns).
    out <- std(data.table::data.table(SNP = "rs1", frequency = 0.2, maf = 0.1))
    expect_equal(names(out), c("SNP", "FRQ", "MAF"))
    expect_equal(out$FRQ, 0.2)

    out <- std(data.table::data.table(SNP = "rs1", FRQ_A = 0.30, FRQ_U = 0.25,
                                      BETA = 0.1))
    expect_equal(names(out), c("SNP", "FRQ_A", "FRQ", "BETA"))
    expect_equal(out$FRQ, 0.25)

    # Neale-lab layout from #138
    out <- std(data.table::data.table(
        variant = "1:15791:C:T", minor_allele = "T", minor_AF = 5.42862e-09,
        AF = 0.99999999457138, beta = -3.56322
    ))
    expect_equal(names(out), c("SNP", "A2", "MINOR_AF", "FRQ", "BETA"))
    expect_equal(anyDuplicated(names(out)), 0)
    expect_equal(out$FRQ, 0.99999999457138)
    expect_equal(out$MINOR_AF, 5.42862e-09)
})

test_that("#45/#48: A0/ALLELE0 becomes A1 and its partner A1/ALLELE1 becomes A2", {
    # A0/A1 means ref/alt, so the input A1 is the effect allele (MSS A2)
    out <- std(data.table::data.table(SNP = "rs1", ALLELE0 = "C", ALLELE1 = "T"))
    expect_equal(names(out), c("SNP", "A1", "A2"))
    expect_equal(c(out$A1, out$A2), c("C", "T"))
    out <- std(data.table::data.table(a0 = "C", a1 = "T"))
    expect_equal(names(out), c("A1", "A2"))
    expect_equal(c(out$A1, out$A2), c("C", "T"))
    # convert_A0 = FALSE leaves the A0 marker so callers can detect it
    out <- std(data.table::data.table(A0 = "C", A1 = "T"), convert_A0 = FALSE)
    expect_equal(names(out), c("A*", "A1"))
})

test_that("#45/#142: a user-supplied mapping_file is applied", {
    # #45: a stand-alone user map (upper-case Uncorrected, as format_sumstats
    # would pass it)
    m <- data.frame(Uncorrected = c("MYSNP", "MYEFF", "MYREF"),
                    Corrected = c("SNP", "A2", "A1"))
    out <- std(data.table::data.table(mySNP = "rs1", myEff = "T", myRef = "C"),
               mapping_file = m)
    expect_equal(names(out), c("SNP", "A2", "A1"))
    expect_equal(c(out$A2, out$A1), c("T", "C"))
    # #142: a row appended to the default map is honoured
    m <- rbind(MungeSumstats:::sumstatsColHeaders,
               data.frame(Uncorrected = "MY_BP_COL", Corrected = "BP"))
    dt <- data.table::data.table(SNP = "rs1", CHR = 1, MY_BP_COL = 100)
    expect_equal(names(std(dt, mapping_file = m)), c("SNP", "CHR", "BP"))
    expect_equal(names(std(dt)), c("SNP", "CHR", "MY_BP_COL"))
})

test_that("#142: REGENIE header (GENPOS, ALLELE0/1, A1FREQ) maps fully", {
    # ~2 s: one untrimmed standardise_header() call, to check that std() is exact
    dt <- data.table::data.table(
        CHROM = 1L, GENPOS = 8490603L, ID = "rs301800", ALLELE0 = "C",
        ALLELE1 = "T", A1FREQ = 0.18, BETA = 0.019, SE = 0.003, P = 1.8e-8
    )
    out <- suppressMessages(MungeSumstats::standardise_header(
        sumstats_dt = data.table::copy(dt), return_list = FALSE
    ))
    expect_equal(names(out),
                 c("CHR", "BP", "SNP", "A1", "A2", "FRQ", "BETA", "SE", "P"))
    # ALLELE1 is REGENIE's effect allele -> A2; A1FREQ is its frequency
    expect_equal(c(out$A1, out$A2), c("C", "T"))
    expect_equal(out$FRQ, 0.18)
    expect_equal(std(dt), out)
})

test_that("#142: in a full REGENIE header A1FREQ, not MAF, becomes FRQ", {
    dt <- data.table::data.table(
        CHROM = 1L, GENPOS = 8490603L, ID = "rs301800", ALLELE0 = "C",
        ALLELE1 = "T", A1FREQ = 0.18, INFO = 1, N = 1000L, TEST = "ADD",
        BETA = 0.019, SE = 0.003, CHISQ = 30, LOG10P = 7.7, P = 1.8e-8,
        MAF = 0.17
    )
    out <- std(dt)
    expect_equal(out$FRQ, 0.18)
    expect_true("MAF" %in% names(out))
})

test_that("#168: allele-specific effect/FRQ names map only via the combinations", {
    # MSS defines BETA/FRQ relative to A2, so A1-specific names are kept out of
    # the base map; infer_effect_column() needs to see them to swap alleles
    m <- MungeSumstats:::sumstatsColHeaders
    expect_false(any(c("BETA1", "A1FRQ", "A1FREQ", "OR_A1", "FRQ_A1") %in%
                         m$Uncorrected))
    cmb <- MungeSumstats:::get_eff_frq_allele_combns()
    hit <- c(BETA1 = "BETA", A1FRQ = "FRQ", A1FREQ = "FRQ", OR_A1 = "OR",
             FREQ.A1.1000G.EUR = "FRQ")
    expect_equal(cmb$CORRECTED[match(names(hit), cmb$UNCORRECTED)],
                 unname(hit))
    out <- std(data.table::data.table(SNP = "rs1", A1FREQ = 0.1, BETA1 = 0.2))
    expect_equal(names(out), c("SNP", "FRQ", "BETA"))
})

test_that("#186: EA/OA map to A2/A1 so A1 is not imputed", {
    input <- data.table::data.table(
        RSID = c("rs1", "rs2"), BETA = c(-0.013, 0.007), SE = c(0.106, 0.10),
        CHISQ = c(0.04, 0.03), P = c(0.77, 0.88), EA = c("A", "T"),
        OA = c("T", "A"), LOG10P = c(0.12, 0.05)
    )
    out <- std(input)
    expect_equal(names(out),
                 c("SNP", "BETA", "SE", "CHISQ", "P", "A2", "A1", "LOG10P"))
    expect_equal(out$A2, c("A", "T"))
    expect_equal(out$A1, c("T", "A"))
})

test_that("#192: risk / non-risk allele headers map to A2 / A1", {
    pairs <- list(
        c("RISK_ALLELE", "NON_RISK_ALLELE"),
        c("risk.allele", "non.risk.allele"),
        c("Risk Allele", "Non Risk Allele"),
        c("RISK-ALLELE", "NON-RISK-ALLELE"),
        c("Risk-Allele", "Non-Risk_Allele"),
        c("RiskAllele", "NON-RISK.ALLELE"),
        c("risk_allele", "Non-Risk Allele")
    )
    for (p in pairs) {
        dt <- data.table::data.table(SNP = c("rs1", "rs2"), x = c("A", "C"),
                                     y = c("G", "T"), BETA = c(0.1, 0.2),
                                     P = 0.5)
        data.table::setnames(dt, c("x", "y"), p)
        out <- std(dt)
        expect_equal(names(out), c("SNP", "A2", "A1", "BETA", "P"),
                     info = paste(p, collapse = " / "))
        expect_equal(out$A2, c("A", "C"))
        expect_equal(out$A1, c("G", "T"))
    }
})

test_that("#99/#195: published header sets map to the exact expected names", {
    sets <- list(
        "#99 Posthuma" = list(
            c("chr", "PosGRCh37", "testedAllele", "otherAllele", "z", "p",
              "N"),
            c("CHR", "BP", "A2", "A1", "Z", "P", "N")),
        "#195 GWAS catalog" = list(
            c("variant_id", "chromosome", "base_pair_location",
              "effect_allele", "other_allele", "odds_ratio", "p_value"),
            c("SNP", "CHR", "BP", "A2", "A1", "OR", "P")),
        "#195 plink" = list(
            c("ID", "CHROM", "POS", "ALT", "REF", "or", "PVAL"),
            c("SNP", "CHR", "BP", "A2", "A1", "OR", "P")),
        "#195 CKDGen" = list(
            c("Chr", "Pos_b37", "RSID", "Allele1", "Allele2", "Freq1",
              "Effect", "StdErr", "P-value", "n_total_sum"),
            c("CHR", "BP", "SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N")),
        # Weight is deliberately unmapped (was -> N before 1.9.18)
        "#195 METAL" = list(
            c("SNP", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq",
              "MaxFreq", "Weight", "Zscore", "P", "Direction", "CHR", "BP"),
            c("SNP", "A1", "A2", "FRQ", "FRQSE", "FRQMIN", "FRQMAX", "WEIGHT",
              "Z", "P", "DIRECTION", "CHR", "BP"))
    )
    for (nm in names(sets)) {
        hdr <- sets[[nm]][[1]]
        dt <- data.table::as.data.table(
            stats::setNames(rep(list(c(1, 2)), length(hdr)), hdr)
        )
        expect_equal(names(std(dt)), sets[[nm]][[2]], info = nm)
    }
})

test_that("#174/#175: format_sumstats normalises a custom mapping_file", {
    # ~9 s: one reference run plus the two custom-map variants (spec needs all 3)
    skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    skip_if_not_installed("BSgenome.Hsapiens.1000genomes.hs37d5")
    run <- function(...) {
        MungeSumstats::format_sumstats(
            path = system.file("extdata", "eduAttainOkbay.txt",
                               package = "MungeSumstats"),
            ref_genome = "GRCh37", dbSNP = 144, on_ref_genome = FALSE,
            strand_ambig_filter = FALSE, bi_allelic_filter = FALSE,
            allele_flip_check = FALSE, return_data = TRUE, ...
        )
    }
    ref <- suppressMessages(run())
    # lower-case Uncorrected entries in a data.frame are upper-cased
    m_low <- MungeSumstats:::sumstatsColHeaders
    m_low$Uncorrected <- tolower(m_low$Uncorrected)
    msgs <- testthat::capture_messages(res_low <- run(mapping_file = m_low))
    expect_true(any(grepl("Non-standard mapping file detected", msgs)))
    expect_equal(res_low, ref)
    # a data.table mapping file passes validate_parameters and is converted
    m_dt <- data.table::as.data.table(MungeSumstats:::sumstatsColHeaders)
    expect_no_error(res_dt <- suppressMessages(run(mapping_file = m_dt)))
    expect_equal(res_dt, ref)
})

test_that("#174/#175: exported standardise_header normalises a custom mapping_file", {
    dt <- data.table::data.table(snp = "rs1", pval = 0.1)
    m_dt <- data.table::as.data.table(MungeSumstats:::sumstatsColHeaders)
    expect_equal(names(std(dt, mapping_file = m_dt)), c("SNP", "P"))
    m_low <- MungeSumstats:::sumstatsColHeaders
    m_low$Uncorrected <- tolower(m_low$Uncorrected)
    expect_equal(names(std(dt, mapping_file = m_low)), c("SNP", "P"))
})
