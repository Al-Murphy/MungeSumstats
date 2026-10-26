# Regression tests for GitHub issues on writing sumstats (write_sumstats, convert_sumstats).
# Issues: #19, #102, #117/#121/#122, #126, #216 (+ the VCF save_path fix), #34 (write_sumstats edge).
# fast stand-in for formatted_example(), whose standardise_header() costs ~2 s
edu_dt <- function() {
    dt <- data.table::fread(
        system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats"),
        nThread = 1
    )
    data.table::setnames(dt, c("SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA",
                               "SE", "P"))
    dt[, CHR := as.character(CHR)]
    suppressMessages(MungeSumstats:::sort_coords(dt))
}

test_that("#102: save_format='openGWAS' writes GWAS-VCF fields", {
    dt <- edu_dt()[1:5]
    dt[, N := 1000L]
    dt[, EXTRA := 1]
    out <- MungeSumstats::write_sumstats(
        sumstats_dt = data.table::copy(dt),
        save_path = tempfile(fileext = ".vcf"),
        write_vcf = TRUE,
        save_format = "openGWAS",
        ref_genome = "GRCh37",
        tabix_index = FALSE,
        return_path = TRUE
    )
    vcf <- VariantAnnotation::readVcf(out)
    # SNP -> RSID in INFO
    expect_equal(names(VariantAnnotation::info(vcf)), "RSID")
    expect_equal(VariantAnnotation::info(vcf)$RSID, dt$SNP)
    # BETA -> ES, P -> LP (-log10), FRQ -> AF, N -> NS; extra columns dropped
    # (asVCF() may add default VRanges FORMAT fields such as AD/DP/FT)
    geno <- VariantAnnotation::geno(vcf)
    expect_true(all(c("NS", "ES", "SE", "LP", "AF") %in% names(geno)))
    expect_false(any(c("EXTRA", "NC", "BETA", "P", "FRQ", "N") %in%
                         names(geno)))
    expect_equal(as.numeric(geno$LP[, 1]), -log10(dt$P), tolerance = 1e-4)
    expect_equal(as.numeric(geno$ES[, 1]), dt$BETA, tolerance = 1e-4)
    expect_equal(as.numeric(geno$SE[, 1]), dt$SE, tolerance = 1e-4)
    expect_equal(as.numeric(geno$AF[, 1]), dt$FRQ, tolerance = 1e-4)
    expect_equal(as.numeric(geno$NS[, 1]), rep(1000, 5))
    expect_error(
        MungeSumstats::write_sumstats(
            sumstats_dt = data.table::copy(dt),
            save_path = tempfile(fileext = ".vcf"), write_vcf = TRUE,
            save_format = "openGWAS", ref_genome = NULL
        ),
        "ref_genome must be provided when write_vcf=TRUE"
    )
    expect_warning(
        MungeSumstats::write_sumstats(
            sumstats_dt = data.table::copy(dt),
            save_path = tempfile(fileext = ".vcf"), write_vcf = TRUE,
            save_format = "openGWAS", ref_genome = "GRCh38"
        ),
        "OpenGWAS is all"
    )
})

test_that("#117/#121/#122: tabix output has integer BP, natural CHR order, no leftovers", {
    skip_if_not_installed("Rsamtools")
    dt <- data.table::data.table(
        SNP = c("rs1", "rs181535630", "rs61573637", "rs6725984", "rs5"),
        CHR = c("10", "2", "2", "2", "1"),
        BP = c(5e7, 71999998, 72000000, 72000148, 1e8),
        A1 = "A", A2 = "G", P = 0.5
    )
    # a double BP is what fwrite prints as 7.2e+07 (tabix then reads "7")
    expect_true(is.double(dt$BP))
    td <- tempfile("t117")
    dir.create(td)
    out <- MungeSumstats::write_sumstats(
        sumstats_dt = data.table::copy(dt),
        save_path = file.path(td, "ss.tsv.gz"),
        tabix_index = TRUE,
        return_path = TRUE
    )
    # #121: only the bgzipped file and its index remain (no residual .tsv)
    expect_equal(out, file.path(td, "ss.tsv.bgz"))
    expect_setequal(list.files(td), c("ss.tsv.bgz", "ss.tsv.bgz.tbi"))
    # #117: no scientific notation, chromosomes in natural order
    txt <- readLines(gzfile(out))
    expect_false(any(grepl("e\\+0", txt)))
    q <- MungeSumstats::read_header(path = out, n = NULL)
    expect_equal(as.character(q$CHR), c("1", "2", "2", "2", "10"))
    expect_equal(as.numeric(q$BP), c(1e8, 71999998, 72000000, 72000148, 5e7))
    # #122: Rsamtools (not seqminer) index answers a region query
    hit <- Rsamtools::scanTabix(
        out, param = GenomicRanges::GRanges("2:72000000-72000000")
    )[[1]]
    expect_length(hit, 1)
    expect_match(hit, "rs61573637")
})

test_that("#126: NA is written as 'NA' in tsv and as an empty field in csv", {
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2"), CHR = c("1", "1"), BP = c(10L, 20L),
        A1 = c("A", "C"), A2 = c("G", "T"), FRQ = c(0.1, NA)
    )
    p <- tempfile(fileext = ".tsv")
    MungeSumstats::write_sumstats(sumstats_dt = data.table::copy(dt),
                                  save_path = p)
    expect_equal(readLines(p)[3], "rs2\t1\t20\tC\tT\tNA")
    p2 <- tempfile(fileext = ".csv")
    MungeSumstats::write_sumstats(sumstats_dt = data.table::copy(dt),
                                  save_path = p2, sep = ",")
    expect_equal(readLines(p2)[3], "rs2,1,20,C,T,")
    skip_if_not_installed("Rsamtools")
    expect_message(
        MungeSumstats::write_sumstats(
            sumstats_dt = data.table::copy(dt),
            save_path = tempfile(fileext = ".tsv.gz"),
            tabix_index = TRUE
        ),
        "Writing uncompressed instead of gzipped to enable tabix indexing"
    )
})

test_that("#19: convert_sumstats returns VRanges (ref=A1, alt=A2), GRanges or the data.table", {
    dt <- edu_dt()
    vr <- MungeSumstats:::convert_sumstats(sumstats_dt = data.table::copy(dt),
                                           return_format = "vranges")
    expect_true(methods::is(vr, "VRanges"))
    expect_length(vr, nrow(dt))
    expect_equal(as.character(VariantAnnotation::ref(vr)), dt$A1)
    expect_equal(as.character(VariantAnnotation::alt(vr)), dt$A2)
    gr <- MungeSumstats:::convert_sumstats(sumstats_dt = data.table::copy(dt),
                                           return_format = "granges")
    expect_true(methods::is(gr, "GRanges"))
    expect_false(methods::is(gr, "VRanges"))
    expect_equal(GenomicRanges::start(gr), dt$BP)
    expect_equal(as.character(GenomicRanges::seqnames(gr)), dt$CHR)
    expect_identical(
        MungeSumstats:::convert_sumstats(sumstats_dt = dt,
                                         return_format = "data.table"),
        dt
    )
})

test_that("#216: .vcf.gz output is BGZF and readable by htslib", {
    p <- MungeSumstats::write_sumstats(
        sumstats_dt = edu_dt(),
        save_path = tempfile(fileext = ".vcf.gz"),
        write_vcf = TRUE,
        return_path = TRUE
    )
    expect_true(file.exists(p))
    con <- file(p, "rb")
    h <- readBin(con, "raw", 18)
    close(con)
    # gzip magic + FEXTRA flag, then the "BC" BGZF subfield
    expect_equal(h[1:4], as.raw(c(0x1f, 0x8b, 0x08, 0x04)))
    expect_equal(rawToChar(h[13:14]), "BC")
    expect_no_error(VariantAnnotation::scanVcfHeader(p))
})

test_that("#216: write_vcf returns an existing path for each VCF suffix and tabix_index", {
    dt <- edu_dt()
    for (ext in c(".vcf", ".vcf.gz", ".vcf.bgz")) {
        for (tbx in c(FALSE, TRUE)) {
            td <- tempfile()
            dir.create(td)
            sp <- file.path(td, paste0("ss", ext))
            p <- suppressMessages(MungeSumstats::write_sumstats(
                sumstats_dt = data.table::copy(dt),
                save_path = sp,
                write_vcf = TRUE,
                tabix_index = tbx,
                return_path = TRUE
            ))
            lab <- sprintf("output for %s, tabix_index=%s", ext, tbx)
            expect_true(file.exists(p), label = lab)
            # only the output (+ index) is left: no temporary uncompressed .vcf
            expect_setequal(list.files(td),
                            basename(c(p, if (tbx) paste0(p, ".tbi"))))
            if (ext == ".vcf.gz" && !tbx) expect_equal(p, sp)
        }
    }
})

test_that("#34: write_sumstats only strips the compression suffix of the file name", {
    skip_if_not_installed("Rsamtools")
    dt <- data.table::data.table(SNP = c("rs1", "rs2"), CHR = "1",
                                 BP = c(10L, 20L), A1 = "A", A2 = "G",
                                 P = 0.5)
    # tabular + tabix: ".gz" also matches "/gz" in the directory part
    td <- tempfile()
    d <- file.path(td, "gz")
    dir.create(d, recursive = TRUE)
    out <- MungeSumstats::write_sumstats(
        sumstats_dt = data.table::copy(dt),
        save_path = file.path(d, "ss.tsv.gz"),
        tabix_index = TRUE,
        return_path = TRUE
    )
    expect_equal(out, file.path(d, "ss.tsv.bgz"))
    expect_false(file.exists(file.path(td, "ss.tsv.bgz")))
    # VCF: a directory name containing ".gz"
    d2 <- file.path(td, "run.gz.d")
    dir.create(d2)
    out2 <- MungeSumstats::write_sumstats(
        sumstats_dt = data.table::copy(dt),
        save_path = file.path(d2, "ss.vcf.gz"),
        write_vcf = TRUE,
        return_path = TRUE
    )
    expect_equal(out2, file.path(d2, "ss.vcf.gz"))
    expect_true(file.exists(out2))
})
