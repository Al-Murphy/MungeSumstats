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

