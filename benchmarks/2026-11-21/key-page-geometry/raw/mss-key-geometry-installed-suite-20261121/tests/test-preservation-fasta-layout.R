test_that("reference FASTA indexes retain exact binary line offsets", {
    # Windows CRLF and a final line without a delimiter must use byte offsets,
    # rather than R's platform-dependent text output or character positions.
    for (ending in c("\n", "\r\n")) {
        for (terminated in c(FALSE, TRUE)) {
            for (width in c(1L, 4L, 7L)) {
                for (short_final in c(FALSE, TRUE)) {
                    d <- tempfile("fasta-layout-")
                    dir.create(d)
                    on.exit(unlink(d, recursive = TRUE), add = TRUE)
                    sequences <- c("ACGTACGTACGTACGTACGTACGTACGT", "GTACGTACGTACGTACGTACGTACGTAC")
                    sequences <- substr(sequences, 1L, width * 3L - as.integer(short_final))
                    fasta <- ""
                    expected <- list()
                    for (i in seq_along(sequences)) {
                        header <- paste0(">", i, " description", ending)
                        offset <- nchar(fasta, type = "bytes") + nchar(header, type = "bytes")
                        starts <- seq.int(1L, nchar(sequences[[i]]), by = width)
                        lines <- substring(sequences[[i]], starts,
                                           pmin(starts + width - 1L, nchar(sequences[[i]])))
                        payload <- paste(lines, collapse = ending)
                        if (i < length(sequences) || terminated) payload <- paste0(payload, ending)
                        fasta <- paste0(fasta, header, payload)
                        expected[[i]] <- c(as.character(i), as.character(nchar(sequences[[i]])),
                                           as.character(offset), as.character(width),
                                           as.character(width + nchar(ending, type = "bytes")))
                    }
                    fasta_path <- file.path(d, "reference.fa")
                    writeBin(charToRaw(fasta), fasta_path)
                    variants <- file.path(d, "variants.tsv")
                    writeBin(charToRaw("chr\tpos\tref\talt\tid\taliases\n"), variants)
                    reference <- prepare_reference(fasta_path, variants, file.path(d, "index"),
                                                   build = "layout-test", version = "1")
                    actual <- read.delim(file.path(reference$path, "fasta.fai"), header = FALSE,
                                         colClasses = "character", check.names = FALSE)
                    expect_identical(unname(as.matrix(actual)), do.call(rbind, expected))
                    # Resolve every sequence position through the public API to
                    # test seeks across line endings and both contig boundaries.
                    input <- do.call(rbind, lapply(seq_along(sequences), function(i) {
                        ref <- strsplit(sequences[[i]], "", fixed = TRUE)[[1L]]
                        alt <- unname(c(A = "C", C = "A", G = "A", T = "C")[ref])
                        data.frame(chr = as.character(i), pos = seq_along(ref),
                                   ref = ref, alt = alt, beta = 0.25)
                    }))
                    result <- harmonise_sumstats(input, reference,
                        columns = list(chr = "chr", pos = "pos", other = "ref",
                                       effect_allele = "alt", effect = "beta",
                                       source_ref = "ref", source_alt = "alt"),
                        source_strand = "+")
                    expect_identical(result$source, input)
                    expect_true(all(result$normalized$status == "sequence_only"))
                    expect_equal(result$normalized$pos, input$pos)
                    expect_identical(result$normalized$ref, input$ref)
                    expect_identical(result$normalized$alt, input$alt)
                }
            }
        }
    }
})

test_that("a single unterminated FASTA sequence line has no delimiter stride", {
    for (ending in c("\n", "\r\n")) {
        d <- tempfile("fasta-single-line-")
        dir.create(d)
        on.exit(unlink(d, recursive = TRUE), add = TRUE)
        fasta <- file.path(d, "reference.fa")
        writeBin(charToRaw(paste0(">1", ending, "ACGT")), fasta)
        variants <- file.path(d, "variants.tsv")
        writeBin(charToRaw("chr\tpos\tref\talt\tid\taliases\n"), variants)
        reference <- prepare_reference(fasta, variants, file.path(d, "index"),
                                       build = "layout-test", version = "1")
        actual <- read.delim(file.path(reference$path, "fasta.fai"), header = FALSE,
                             colClasses = "character")
        expect_identical(unname(unlist(actual)),
                         c("1", "4", as.character(2L + nchar(ending, type = "bytes")), "4", "4"))
    }
})
