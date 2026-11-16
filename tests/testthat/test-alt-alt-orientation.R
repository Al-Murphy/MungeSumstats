# A1/A2 are measured alleles, not a declaration of source REF or strand.
# A single REF-ALT match does not resolve a competing, coherent ALT-ALT pair.
.alt_orientation_reference <- function(rows, sequences = c(`1` = "AAAAAA", `2` = "AAAAAA")) {
    root <- tempfile("alt-orientation-")
    dir.create(root)
    fasta <- file.path(root, "reference.fa")
    writeLines(unlist(Map(function(name, sequence) c(paste0(">", name), sequence),
                          names(sequences), unname(sequences))), fasta)
    variants <- file.path(root, "variants.tsv")
    writeLines(c("chr\tpos\tref\talt\tid\taliases", rows), variants)
    list(root = root, reference = prepare_reference(
        fasta, variants, file.path(root, "index"), build = "test", version = "alt-orientation"))
}

.alt_orientation_query <- function(other = "A", effect = "G", chr = "1", pos = 1,
                                   id = "shared") {
    data.frame(CHR = chr, BP = pos, SNP = id, A1 = other, A2 = effect,
               BETA = "0.300", FRQ = "0.200", P = "1e-1000",
               untouched = "original-token", stringsAsFactors = FALSE)
}

.alt_orientation_expect_unresolved <- function(normalized, status = "ambiguous") {
    expect_identical(normalized$status, rep(status, nrow(normalized)))
    expect_true(all(is.na(normalized[c("chr", "pos", "ref", "alt", "id", "effect", "eaf")])))
    expect_true(all(is.na(normalized$swapped)))
    expect_true(all(is.na(normalized$complemented)))
}

test_that("quad-allelic source contrasts remain ambiguous across all orientations", {
    fixture <- .alt_orientation_reference("1\t1\tA\tC,G,T\tshared\told_shared")
    on.exit(unlink(fixture$root, recursive = TRUE))
    # A/G and A/C each have a reverse-strand ALT-ALT interpretation; their
    # swaps and reverse complements carry the same unresolved alternatives.
    input <- .alt_orientation_query(c("A", "G", "T", "C", "A", "C", "T", "G"),
                                    c("G", "A", "C", "T", "C", "A", "G", "T"))
    actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
    expect_identical(actual$source, input)
    expect_equal(nrow(actual$normalized), nrow(input))
    expect_equal(actual$normalized$candidate_count, rep(3, nrow(input)))
    .alt_orientation_expect_unresolved(actual$normalized)
    input$SNP <- "old_shared"
    alias <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
    expect_identical(alias$source, input)
    fields <- setdiff(names(actual$normalized), "entity_id")
    expect_identical(alias$normalized[fields], actual$normalized[fields])
    input$SNP <- NA_character_
    expect_identical(harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")$normalized[fields],
                     actual$normalized[fields])
})

test_that("trusted source REF ALT and strand still resolve a quad-allelic event", {
    fixture <- .alt_orientation_reference("1\t1\tA\tC,G,T\tshared\t.")
    on.exit(unlink(fixture$root, recursive = TRUE))
    columns <- list(chr = "CHR", pos = "BP", id = "SNP", other = "A1",
                    effect_allele = "A2", effect = "BETA", eaf = "FRQ",
                    source_ref = "REF", source_alt = "ALT", source_strand = "STRAND")
    input <- .alt_orientation_query(c("A", "G", "T", "C"), c("G", "A", "C", "T"))
    input$REF <- c("A", "A", "T", "T")
    input$ALT <- c("G", "G", "C", "C")
    input$STRAND <- c("+", "+", "-", "-")
    actual <- harmonise_sumstats(input, fixture$reference, columns = columns, frequency_type = "EAF")
    expect_identical(actual$source, input)
    expect_identical(actual$normalized$status, c("match", "swap", "complement", "swapped_complement"))
    expect_identical(actual$normalized$ref, rep("A", 4))
    expect_identical(actual$normalized$alt, rep("G", 4))
    expect_equal(actual$normalized$effect, c(.3, -.3, .3, -.3))
    expect_equal(actual$normalized$eaf[c(1, 3)], c(.2, .2))
    expect_true(all(is.na(actual$normalized$eaf[c(2, 4)])))
    input$STRAND <- "unknown"
    unknown <- harmonise_sumstats(input, fixture$reference, columns = columns)$normalized
    .alt_orientation_expect_unresolved(unknown, "source_strand_unknown")
})

test_that("split biallelic source records retain coherent orientation ambiguity regardless of order", {
    rows <- c("1\t1\tA\tC\tshared\t.", "1\t1\tA\tG\tshared\t.",
              "1\t1\tA\tT\tshared\t.")
    orders <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
    input <- .alt_orientation_query(c("A", "G", "T", "C"), c("G", "A", "C", "T"))
    previous <- NULL
    for (order in orders) {
        fixture <- .alt_orientation_reference(rows[order])
        on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
        actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF", batch_size = 1)
        expect_identical(actual$source, input)
        .alt_orientation_expect_unresolved(actual$normalized)
        if (!is.null(previous)) expect_identical(actual$normalized, previous)
        previous <- actual$normalized
    }
})

test_that("indel and multi-base ALT pairs can compete with a single canonical event", {
    indel <- .alt_orientation_reference("1\t1\tA\tAT,T\tshared\t.")
    on.exit(unlink(indel$root, recursive = TRUE), add = TRUE)
    indel_input <- .alt_orientation_query(c("A", "AT", "T", "AT"), c("AT", "A", "AT", "T"))
    actual <- harmonise_sumstats(indel_input, indel$reference, frequency_type = "EAF")
    expect_identical(actual$source, indel_input)
    .alt_orientation_expect_unresolved(actual$normalized)

    mnv <- .alt_orientation_reference("1\t1\tAC\tAT,GT\tshared\t.", c(`1` = "ACGAAA"))
    on.exit(unlink(mnv$root, recursive = TRUE), add = TRUE)
    mnv_input <- .alt_orientation_query(c("AC", "AT", "GT", "AT"), c("AT", "AC", "AT", "GT"))
    actual <- harmonise_sumstats(mnv_input, mnv$reference, frequency_type = "EAF")
    expect_identical(actual$source, mnv_input)
    .alt_orientation_expect_unresolved(actual$normalized)
})

test_that("biallelic and noncompeting triallelic directions retain resolved effects", {
    fixtures <- list(
        .alt_orientation_reference("1\t1\tA\tG\tshared\t."),
        .alt_orientation_reference("1\t1\tA\tC,G\tshared\t."))
    input <- .alt_orientation_query(c("A", "G", "T", "C"), c("G", "A", "C", "T"))
    for (i in seq_along(fixtures)) {
        fixture <- fixtures[[i]]
        on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
        actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
        expect_identical(actual$source, input)
        expect_identical(actual$normalized$status, c("match", "swap", "complement", "swapped_complement"))
        expect_equal(actual$normalized$effect, c(.3, -.3, .3, -.3))
        expect_equal(actual$normalized$eaf[c(1, 3)], c(.2, .2))
        if (i == 1) expect_equal(actual$normalized$eaf[c(2, 4)], c(.8, .8))
        else expect_true(all(is.na(actual$normalized$eaf[c(2, 4)])))
    }
    alt_pair <- harmonise_sumstats(.alt_orientation_query("C", "G"), fixtures[[2]]$reference)$normalized
    .alt_orientation_expect_unresolved(alt_pair, "alt_alt_contrast")
    pal <- .alt_orientation_reference("1\t1\tA\tT\tshared\t.")
    on.exit(unlink(pal$root, recursive = TRUE), add = TRUE)
    .alt_orientation_expect_unresolved(harmonise_sumstats(.alt_orientation_query("A", "T"), pal$reference)$normalized)
})

test_that("incoherent or invalid ALT contexts cannot obscure a unique canonical hit", {
    controls <- list(
        # T and C must not be pooled across different original REF spans.
        c("1\t1\tA\tG,T\tshared\t.", "1\t1\tAA\tC\tother\t."),
        # A genuine T/C pair at another position is not the queried site.
        c("1\t1\tA\tG\tshared\t.", "1\t3\tA\tT,C\tshared\t."),
        # Original site records survive preparation even when REF is unsupported
        # or inconsistent with FASTA. Such records cannot validate an ALT pair.
        c("1\t1\tA\tG\tshared\t.", "1\t1\tN\tT,C\tother\t."),
        c("1\t1\tA\tG\tshared\t.", "1\t1\tG\tT,C\tother\t."))
    input <- .alt_orientation_query()
    for (rows in controls) {
        fixture <- .alt_orientation_reference(rows)
        on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
        actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
        expect_identical(actual$source, input)
        expect_identical(actual$normalized$status, "match")
        expect_identical(actual$normalized$ref, "A")
        expect_identical(actual$normalized$alt, "G")
        expect_equal(actual$normalized$effect, .3)
        expect_equal(actual$normalized$eaf, .2)
    }
    foreign <- .alt_orientation_reference(c("1\t1\tA\tG\tshared\t.", "2\t1\tA\tC,G,T\tshared\t."))
    on.exit(unlink(foreign$root, recursive = TRUE), add = TRUE)
    input$BP <- NA_real_
    actual <- harmonise_sumstats(input, foreign$reference, frequency_type = "EAF")
    expect_identical(actual$normalized$status, "match")
    expect_identical(actual$normalized$chr, "1")
    expect_equal(actual$normalized$effect, .3)
})

test_that("a bounded site-context lookup cannot silently select a unique hit on overflow", {
    fixture <- .alt_orientation_reference(c("1\t1\tA\tG\tshared\t.",
                                            "1\t1\tN\tC,T\tother\t.",
                                            "1\t1\tNN\tC,T\tother\t."))
    on.exit(unlink(fixture$root, recursive = TRUE))
    expect_error(harmonise_sumstats(.alt_orientation_query(), fixture$reference,
                                   frequency_type = "EAF", candidate_limit = 1),
                 "candidate limit exceeded")
})

test_that("streaming preserves original bytes and every ambiguous source row", {
    fixture <- .alt_orientation_reference("1\t1\tA\tC,G,T\tshared\t.")
    on.exit(unlink(fixture$root, recursive = TRUE))
    input <- .alt_orientation_query(c("A", "G", "T", "C"), c("G", "A", "C", "T"))
    path <- file.path(fixture$root, "input.tsv")
    write.table(input, path, sep = "\t", row.names = FALSE, quote = FALSE)
    original_bytes <- readBin(path, "raw", n = file.info(path)$size)
    parsed <- read.delim(path, colClasses = "character", check.names = FALSE)
    resident <- harmonise_sumstats(parsed, fixture$reference, frequency_type = "EAF")
    for (batch_size in c(1, 3, 100)) {
        output <- harmonise_file(path, fixture$reference,
                                 file.path(fixture$root, paste0("edition-", batch_size)),
                                 frequency_type = "EAF", batch_size = batch_size)
        expect_equal(output$source_rows, nrow(input))
        expect_equal(output$normalized_rows, nrow(input))
        expect_identical(readBin(file.path(output$path, "source-artifact"), "raw", n = length(original_bytes)),
                         original_bytes)
        expect_identical(readLines(file.path(output$path, "source.tsv")), readLines(path))
        normalized <- read.delim(file.path(output$path, "normalized.tsv"), check.names = FALSE)
        expect_identical(normalized$status, resident$normalized$status)
        expect_identical(normalized$status, rep("ambiguous", nrow(input)))
        expect_true(all(is.na(normalized[c("chr", "pos", "ref", "alt", "id", "effect", "eaf")])))
        expect_equal(normalized$source_row, seq_len(nrow(input)))
    }
    expect_identical(resident$source, parsed)
})
