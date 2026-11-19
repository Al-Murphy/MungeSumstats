legacy_chr_style <- function(chr, style) {
    GenomeInfoDb::mapSeqlevels(as.character(chr), style = style)
}

test_that("chromosome style mapping preserves values, row order and names", {
    chromosomes <- c(as.character(1:22), "X", "Y", "MT")
    fixtures <- list(
        rep("22", 128L),
        c(rev(chromosomes), chromosomes[c(22, 1, 25, 23, 2, 24, 22)]),
        c("chr22", "chrX", "chrM", "chrY", "chr22", "chr1"),
        c("ch22", "chX", "chMT", "chY", "ch22", "ch1"),
        factor(c("MT", "22", "X", "22")),
        stats::setNames(c("2", "1", "2", "MT"), letters[1:4]),
        chromosomes
    )
    for (style in c("NCBI", "Ensembl", "UCSC", "dbSNP")) {
        for (chr in fixtures) {
            expect_identical(map_chr_style(chr, style),
                             legacy_chr_style(chr, style))
        }
    }

    chr <- c("22", "1", "22")
    expect_identical(map_chr_style(chr, "UCSC"),
                     stats::setNames(c("chr22", "chr1", "chr22"), chr))
})

test_that("best chromosome map depends on the frequency of competing styles", {
    ncbi_majority <- c("1", "chr1", "1", "1")
    ucsc_majority <- c("1", "chr1", "chr1", "chr1")
    for (chr in list(ncbi_majority, ucsc_majority)) {
        expected <- legacy_chr_style(chr, "dbSNP")
        expect_identical(map_chr_style(chr, "dbSNP"), expected)
    }
    expect_identical(unname(legacy_chr_style(ncbi_majority, "dbSNP")),
                     c("ch1", NA, "ch1", "ch1"))
    expect_identical(unname(legacy_chr_style(ucsc_majority, "dbSNP")),
                     c(NA, "ch1", "ch1", "ch1"))

    other_majority <- c("1", "2", rep("chr3", 4L))
    expected <- legacy_chr_style(other_majority, "dbSNP")
    expect_identical(map_chr_style(other_majority, "dbSNP"), expected)
    contigs <- unique(other_majority)
    unweighted <- legacy_chr_style(contigs, "dbSNP")
    expect_false(identical(unweighted[match(other_majority, contigs)], expected))

    tie <- c("chr1", "1", "chr1", "1")
    expected_tie <- legacy_chr_style(tie, "dbSNP")
    expect_identical(dim(expected_tie), c(2L, 4L))
    expect_identical(map_chr_style(tie, "dbSNP"), expected_tie)
})

test_that("chromosome mapping preserves unmapped contigs and zero-map matrices", {
    contigs <- c("GL000192.1", "chr6_GL000250v2_alt",
                 "chr22_KI270731v1_random", "unknown_contig")
    for (chr in list(c("22", contigs, "22"),
                     c("22", NA_character_, "", "22", ""),
                     c(contigs, contigs))) {
        for (style in c("NCBI", "Ensembl", "UCSC", "dbSNP")) {
            expect_identical(map_chr_style(chr, style),
                             legacy_chr_style(chr, style))
        }
    }
})

test_that("chromosome mapping preserves legacy warnings and errors", {
    mapping_conditions <- function(mapper, chr, style) {
        warnings <- list()
        result <- withCallingHandlers(
            tryCatch(mapper(chr, style), error = function(error) {
                list(error = conditionMessage(error), class = class(error))
            }),
            warning = function(warning) {
                warnings[[length(warnings) + 1L]] <<-
                    list(message = conditionMessage(warning), class = class(warning))
                invokeRestart("muffleWarning")
            }
        )
        list(result = result, warnings = warnings)
    }
    expect_identical(mapping_conditions(map_chr_style, character(), "UCSC"),
                     mapping_conditions(legacy_chr_style, character(), "UCSC"))

    chr <- c("22", "MT", "22")
    for (style in list("unsupported-chromosome-style", c("UCSC", "NCBI"),
                       NA_character_, "", NULL)) {
        expect_identical(mapping_conditions(map_chr_style, chr, style),
                         mapping_conditions(legacy_chr_style, chr, style))
    }
})
