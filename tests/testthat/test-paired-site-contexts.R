.paired_site_reference <- function(rows) {
    root <- tempfile("paired-site-")
    dir.create(root)
    fasta <- file.path(root, "reference.fa")
    variants <- file.path(root, "variants.tsv")
    writeLines(c(">1", "CAAAC", ">2", "GGGGG", ">3", "ATGCG"), fasta)
    writeLines(c("chr\tpos\tref\talt\tid\taliases", rows), variants)
    list(root = root, reference = prepare_reference(
        fasta, variants, file.path(root, "index"), build = "test", version = "paired-sites"))
}

.paired_site_query <- function(chr = "1", pos = 4, id = "shared", other = "AA", effect = "AAA") {
    data.frame(CHR = chr, BP = pos, SNP = id, A1 = other, A2 = effect,
               BETA = 0.5, FRQ = 0.3, P = "1e-1000", stringsAsFactors = FALSE)
}

test_that("separate original sites cannot invent an ALT-ALT diagnostic through record order", {
    rows <- c("1\t1\tC\tT\tpad\t.", "1\t2\tAA\tAAA\tshared\t.",
              "1\t4\tA\tAA\tshared\t.")
    orders <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
    input <- .paired_site_query()
    for (order in orders) {
        fixture <- .paired_site_reference(rows[order])
        on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
        actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
        expect_identical(actual$source, input)
        expect_identical(actual$normalized$status, "allele_mismatch", info = paste(order, collapse = ","))
        expect_true(all(is.na(actual$normalized[c("pos", "ref", "alt", "id", "effect", "eaf")])))
    }
})

test_that("an ALT pair in one original REF span remains a genuine ALT-ALT diagnostic", {
    fixture <- .paired_site_reference("2\t3\tG\tA,C\tmulti\told_multi")
    on.exit(unlink(fixture$root, recursive = TRUE))
    input <- .paired_site_query("2", 3, "old_multi", "A", "C")
    actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
    expect_identical(actual$normalized$status, "alt_alt_contrast")
    expect_identical(actual$source, input)
    expect_true(is.na(actual$normalized$effect))
    expect_true(is.na(actual$normalized$eaf))
})

test_that("different REF spans at one position cannot pool ALT strings", {
    fixture <- .paired_site_reference(c("2\t3\tG\tA\tspan\told_span",
                                       "2\t3\tGG\tCC\tspan\told_span"))
    on.exit(unlink(fixture$root, recursive = TRUE))
    actual <- harmonise_sumstats(.paired_site_query("2", 3, "old_span", "A", "CC"),
                                fixture$reference, frequency_type = "EAF")
    expect_identical(actual$normalized$status, "allele_mismatch")
    expect_true(all(is.na(actual$normalized[c("pos", "effect", "eaf")])))
})

test_that("an ALT pair at a different original position cannot classify the input site", {
    fixture <- .paired_site_reference(c("1\t2\tAA\tAAA,AAAA\tshared\t.",
                                       "1\t4\tA\tAA\tshared\t."))
    on.exit(unlink(fixture$root, recursive = TRUE))
    actual <- harmonise_sumstats(.paired_site_query(other = "AAA", effect = "AAAA"),
                                fixture$reference, frequency_type = "EAF")
    expect_identical(actual$normalized$status, "allele_mismatch")
    expect_true(all(is.na(actual$normalized[c("pos", "effect", "eaf")])))
})

test_that("an ID-only lookup cannot borrow an ALT pair from a different known chromosome", {
    fixture <- .paired_site_reference(c("2\t3\tG\tA,T\tshared\t.",
                                       "3\t3\tG\tC\tshared\t."))
    on.exit(unlink(fixture$root, recursive = TRUE))
    actual <- harmonise_sumstats(.paired_site_query("3", NA_real_, "shared", "A", "T"),
                                fixture$reference, frequency_type = "EAF")
    expect_identical(actual$normalized$status, "allele_mismatch")
    expect_true(all(is.na(actual$normalized[c("pos", "effect", "eaf")])))
})

test_that("an unsupported original REF is context rather than validated ALT-ALT evidence", {
    fixture <- .paired_site_reference(c("2\t3\tG\tC\tvalid\t.",
                                       "2\t3\tN\tA,T\t.\t."))
    on.exit(unlink(fixture$root, recursive = TRUE))
    actual <- harmonise_sumstats(.paired_site_query("2", 3, "valid", "A", "T"),
                                fixture$reference, frequency_type = "EAF")
    expect_identical(actual$normalized$status, "allele_mismatch")
    expect_true(all(is.na(actual$normalized[c("pos", "effect", "eaf")])))
})

test_that("symbolic context preserves real ALT pairs and conservative resolved frequency", {
    fixture <- .paired_site_reference("2\t3\tG\tA,T,<DEL>\tmulti\told_multi")
    on.exit(unlink(fixture$root, recursive = TRUE))
    contrast <- harmonise_sumstats(.paired_site_query("2", 3, "old_multi", "A", "T"),
                                  fixture$reference, frequency_type = "EAF")
    expect_identical(contrast$normalized$status, "alt_alt_contrast")
    swap <- harmonise_sumstats(.paired_site_query("2", 3, "old_multi", "A", "G"),
                              fixture$reference, frequency_type = "EAF")
    expect_identical(swap$normalized$status, "swap")
    expect_identical(swap$normalized$effect, -0.5)
    expect_true(is.na(swap$normalized$eaf))
})

test_that("biallelic direction and frequency remain resolved across all four orientations", {
    fixture <- .paired_site_reference("2\t3\tG\tA\tbi\told_bi")
    on.exit(unlink(fixture$root, recursive = TRUE))
    input <- data.frame(CHR = "2", BP = 3, SNP = "old_bi", A1 = c("G", "A", "C", "T"),
                        A2 = c("A", "G", "T", "C"), BETA = 0.5, FRQ = 0.3)
    actual <- harmonise_sumstats(input, fixture$reference, frequency_type = "EAF")
    expect_identical(actual$source, input)
    expect_identical(actual$normalized$status, c("match", "swap", "complement", "swapped_complement"))
    expect_identical(actual$normalized$effect, c(0.5, -0.5, 0.5, -0.5))
    expect_equal(actual$normalized$eaf, c(0.3, 0.7, 0.3, 0.7))
})
