.stream_api <- function(name) get(name, envir = environment(harmonise_file))

test_that("file publication no longer materializes legacy row lists", {
    d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
    fasta <- file.path(d, "fasta"); sites <- file.path(d, "sites")
    writeLines(c(">1", "ACGT"), fasta)
    writeLines(c("chr\tpos\tref\talt\tid\taliases", "1\t2\tC\tA\tr\t."), sites)
    ref <- prepare_reference(fasta, sites, file.path(d, "ref"), "synthetic", "toy")
    input <- file.path(d, "in.tsv")
    writeLines(c("CHR\tBP\tA1\tA2\tSNP", "1\t2\tC\tA\tr", "1\t2\tA\tC\tr"), input)
    run <- harmonise_file
    env <- new.env(parent = environment(run))
    env$mss_parser_feed_cpp <- function(...) stop("legacy row-list parser used for file publication")
    environment(run) <- env
    result <- run(input, ref, file.path(d, "edition"), batch_size = 2L, batch_bytes = 8L)
    expect_equal(result$source_rows, 2)
    normalized <- read.delim(file.path(result$path, "normalized.tsv"))
    expect_identical(normalized$status, c("match", "swap"))
})

# Drive exactly the cursor/EOF protocol used on a real chunked connection.
.stream_batches <- function(chunks, width, keep, rows, budget, cap = 10000L) {
    batch <- .stream_api("mss_parser_batch_cpp")
    parser <- .stream_api("mss_parser_open_cpp")(",", cap)
    header <- NULL
    batches <- list()
    for (chunk in c(chunks, list(raw()))) {
        offset <- 0L
        eof <- !length(chunk)
        # A batch can be returned without consuming a byte (a deferred record).
        for (attempt in seq_len(4L + length(chunk) + rows)) {
            reading_header <- is.null(header)
            value <- batch(parser, chunk, offset, eof,
                           if (reading_header) 1L else rows, budget,
                           if (reading_header) 0L else width,
                           if (reading_header) integer() else keep)
            stopifnot(identical(names(value), c("state", "offset", "nrow", "columns")),
                      value$offset >= offset, value$offset <= length(chunk))
            offset <- value$offset
            if (value$state != "batch") {
                stopifnot(value$nrow == 0L, identical(value$columns, list()))
                if (value$state == "done") {
                    stopifnot(eof)
                    return(list(header = header, batches = batches))
                }
                stopifnot(identical(value$state, "need_input"), !eof, offset == length(chunk))
                break
            }
            if (reading_header) {
                stopifnot(value$nrow == 1L, length(value$columns) == 1L)
                header <- value$columns[[1L]]
            } else {
                stopifnot(value$nrow >= 1L, value$nrow <= rows,
                          length(value$columns) == length(keep), is.null(names(value$columns)),
                          all(lengths(value$columns) == value$nrow))
                batches[[length(batches) + 1L]] <- list(nrow = value$nrow, columns = value$columns)
            }
            if (!eof && offset == length(chunk)) break
        }
    }
    stop("Parser failed to drain EOF")
}

.stream_groups <- function(records, keep, rows, budget) {
    answer <- list()
    next_row <- 1L
    while (next_row <= length(records)) {
        first <- next_row
        used <- 0
        while (next_row <= length(records) && next_row - first < rows) {
            size <- sum(nchar(records[[next_row]], type = "bytes")) + length(records[[next_row]])
            if (next_row > first && used + size > budget) break
            used <- used + size
            next_row <- next_row + 1L
        }
        selected <- records[first:(next_row - 1L)]
        answer[[length(answer) + 1L]] <- list(nrow = length(selected),
            columns = lapply(keep, function(j) vapply(selected, `[[`, character(1), j)))
    }
    answer
}

test_that("column batches match row decoding and row/byte membership at every split", {
    tokens <- list(c("A", "B", "C"), c("1", "x", ""), c("", "", "a"),
                   c('"', "雪\r\n,β", "three"), c("4", strrep("x", 64), "4"),
                   c("5", "x", ""), c("6", "x", ""), c("7", "x", ""))
    encode <- function(row) paste0('"', gsub('"', '""', row, fixed = TRUE), '"', collapse = ",")
    bytes <- charToRaw(enc2utf8(paste(vapply(tokens, encode, character(1)), collapse = "\r\n")))
    legacy <- .stream_api("mss_parser_feed_cpp")(.stream_api("mss_parser_open_cpp")(",", 10000L),
                                                  bytes, TRUE)
    expect_identical(legacy, tokens)
    check <- function(chunks, keep, rows, budget) {
        got <- .stream_batches(chunks, 3L, keep, rows, budget)
        expect_identical(got$header, legacy[[1]])
        expect_identical(got$batches, .stream_groups(legacy[-1L], keep, rows, budget))
    }
    for (i in seq_len(length(bytes) - 1L)) {
        check(list(bytes[seq_len(i)], bytes[-seq_len(i)]), c(1L, 3L), 3L, 10L)
    }
    for (keep in list(integer(), 2L, 1:3)) {
        for (budget in c(1L, 5L, 9L, 10L, 1000L)) {
            check(as.list(bytes), keep, 2L, budget)
        }
    }
    expect_null(.stream_batches(list(), 2L, 1L, 2L, 4L)$header)
    header_only <- .stream_batches(list(charToRaw("A,B")), 2L, 1L, 2L, 1L)
    expect_identical(header_only, list(header = c("A", "B"), batches = list()))
})

test_that("collectors reject mixed use and invalid source records even in discarded columns", {
    open <- .stream_api("mss_parser_open_cpp")
    feed <- .stream_api("mss_parser_feed_cpp")
    batch <- .stream_api("mss_parser_batch_cpp")
    parser <- open(",", 100L)
    feed(parser, charToRaw("A"), FALSE)
    expect_error(batch(parser, raw(), 0L, TRUE, 1L, 8L, 0L, integer()))
    parser <- open(",", 100L)
    batch(parser, charToRaw("A"), 0L, FALSE, 1L, 8L, 0L, integer())
    expect_error(feed(parser, raw(), TRUE))
    expect_error(batch(open(",", 100L), charToRaw("A"), 0L, TRUE, 1L, 8L, 0L, integer()))
    for (invalid in list(c(1L, 1L), c(2L, 1L), 3L, 0L, NA_integer_)) {
        expect_error(.stream_batches(list(charToRaw("A,B\n1,2")), 2L, invalid, 1L, 8L))
    }
    cases <- list(list(charToRaw("A,B\n1,x,y"), "fields"),
                  list(c(charToRaw("A,B\n1,x"), as.raw(0L)), "NUL"),
                  list(charToRaw('A,B\n1,"x"x'), "quote"),
                  list(charToRaw('A,B\n1,"x'), "quoted"),
                  list(charToRaw('A,B\n1,q"x'), "quote"),
                  list(charToRaw("A,B\n1,1234567"), "max_record_bytes"))
    for (case in cases) {
        expect_error(.stream_batches(as.list(case[[1L]]), 2L, 1L, 3L, 9L, cap = 8L), case[[2L]])
    }
})

.stream_baseline_write <- function(table, con, header) {
    encoded <- lapply(table, function(value) {
        if (!is.numeric(value)) return(as.character(value))
        text <- sprintf("%.17g", value)
        text[is.na(value)] <- NA_character_
        text
    })
    utils::write.table(as.data.frame(encoded, stringsAsFactors = FALSE, check.names = FALSE),
        con, sep = "\t", quote = TRUE, row.names = FALSE, col.names = header,
        na = "NA", qmethod = "double")
}

.stream_output <- function(writer, pieces, ...) {
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    con <- file(path, open = "wt")
    tryCatch({
        for (i in seq_along(pieces)) writer(pieces[[i]], con, i == 1L, ...)
    }, finally = close(con))
    readBin(path, "raw", n = file.info(path)$size)
}

test_that("bounded writer matches exact baseline bytes across block boundaries", {
    encode <- .stream_api("mss_write_normalized_cpp")
    numbers <- c(-0.0, 0.0, 1, 1 + .Machine$double.eps, 1 - .Machine$double.eps / 2,
                 .Machine$double.xmin * .Machine$double.eps, .Machine$double.xmax,
                 1e-4, 1e-5, 1e16, 1e17, 2^31, 2^53 - 1, Inf, -Inf, NaN, NA_real_)
    tab <- data.frame(pos = numbers, effect = rev(numbers),
                      status = rep_len(c("", "NA", NA_character_, '"\t\r\n', "β雪", 'x"y'), length(numbers)),
                      swapped = rep_len(c(TRUE, FALSE, NA), length(numbers)),
                      candidate_count = rep_len(c(NA_integer_, 0L, 23L, .Machine$integer.max), length(numbers)))
    # A single field larger than many output blocks, with escape and UTF-8 boundaries.
    tab$status[[2L]] <- strrep('雪"\t\r\n', 80L)
    for (pieces in list(list(tab[FALSE, ]), list(tab), list(tab[1:4, ], tab[5:17, ]))) {
        expected <- .stream_output(.stream_baseline_write, pieces)
        for (block in c(1L, 7L, 257L, 262144L)) {
            expect_identical(.stream_output(encode, pieces, block_bytes = block), expected)
        }
    }
    opts <- options(encoding = "UTF-8")
    on.exit(options(opts), add = TRUE)
    expect_identical(.stream_output(encode, list(tab), block_bytes = 1L),
                     .stream_output(.stream_baseline_write, list(tab)))
})

test_that("file and resident batches retain source contracts, classes and edition bytes", {
    d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
    fasta <- file.path(d, "fasta"); sites <- file.path(d, "sites")
    writeLines(c(">1", "CAAAAACGTACGT"), fasta)
    writeLines(c("chr\tpos\tref\talt\tid\taliases", "1\t8\tG\tA,T\trs8\tmerged8",
                 "1\t4\tAA\tA\trsDel\t.", "1\t10\tA\tC\trs10\t."), sites)
    x <- data.frame(chr = factor("chr1"), pos = factor(c(8, 10, 4, 8)), id = c("merged8", "rs10", "rsDel", "merged8"),
                    other = c("T", "A", "AA", "G"), measured = c("G", "C", "A", "T"),
                    beta = c(1 + .Machine$double.eps, -.125, .25, 3), eaf = c(.2, .25, .5, .75),
                    ref = c("G", "A", "AA", "G"), alt = c("T", "C", "A", "T"),
                    kind = c("", "", "", "gene_burden"),
                    name = as.Date("2021-01-01") + 0:3,
                    spare = strrep("x", 600), stringsAsFactors = FALSE)
    x <- x[c(4L, 2L, 1L, 3L, 1L, 2L), ]
    cols <- list(chr = "chr", pos = "pos", id = "id", other = "other", effect_allele = "measured",
                 effect = "beta", eaf = "eaf", source_ref = "ref", source_alt = "alt",
                 entity_type = "kind", entity_id = "name")
    chars <- as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
    input <- file.path(d, "input.csv")
    write.table(chars, input, sep = ",", quote = TRUE, qmethod = "double", row.names = FALSE)
    original <- readBin(input, "raw", n = file.info(input)$size)
    for (fmt in 2:3) {
        ref <- prepare_reference(fasta, sites, file.path(d, paste0("r", fmt)), "syn", "toy",
                  contig_aliases = c(chr1 = "1"), reference_format = paste0("mss-preservation-", fmt))
        args <- list(reference = ref, columns = cols, source_strand = "+", frequency_type = "EAF")
        resident <- do.call(harmonise_sumstats, c(list(sumstats = x, batch_size = 2L), args))
        expect_identical(resident$source, x)
        expected <- resident$normalized
        expect_identical(expected$status, c("specialist_preserved", "match", "swap", "match", "swap", "match"))
        expect_identical(expected$entity_id, as.character(x$name))
        expect_identical(expected$effect, c(NA_real_, -.125, -1, .25, -1, -.125))
        expect_identical(expected$eaf, c(NA_real_, .25, NA_real_, .5, NA_real_, .25))
        expect_identical(expected, do.call(harmonise_sumstats, c(list(sumstats = chars), args))$normalized)
        out <- do.call(harmonise_file, c(list(input = input, path = file.path(d, paste0("out", fmt)),
            delimiter = ",", batch_size = 3L, batch_bytes = 32L,
            acquisition_receipt = list(rows = nrow(x), bytes = length(original), md5 = unname(tools::md5sum(input)))), args))
        for (file in c("source.csv", "source-artifact")) {
            expect_identical(readBin(file.path(out$path, file), "raw", length(original) + 1L), original)
        }
        normalized <- file.path(out$path, "normalized.tsv")
        expect_identical(readBin(normalized, "raw", file.info(normalized)$size),
                         .stream_output(.stream_baseline_write, list(expected)))
        expect_identical(unname(tools::md5sum(file.path(out$path, names(out$checksums)))), unname(out$checksums))
    }
    failed <- file.path(d, "failed")
    writeLines(c("A,A", "1,x"), input)
    expect_error(harmonise_file(input, ref, failed, delimiter = ","), "Header")
    expect_false(dir.exists(failed))
    writeBin(raw(), input)
    expect_error(harmonise_file(input, ref, failed, delimiter = ","), "header")
    expect_false(dir.exists(failed))
    writeBin(charToRaw("other,unused\r\n1,x\r\n2,y"), input)
    expect_error(harmonise_file(input, ref, failed, delimiter = ",", acquisition_receipt = list(rows = 9L)), "receipt")
    expect_false(dir.exists(failed))
    absent <- harmonise_file(input, ref, file.path(d, "absent"), delimiter = ",", columns = list(chr = "missing"))
    expect_equal(absent$source_rows, 2)
    con <- gzfile(input, open = "wb"); writeBin(original, con); close(con)
    compressed <- readBin(input, "raw", file.info(input)$size)
    writeBin(head(compressed, -8L), input)
    expect_error(harmonise_file(input, ref, failed, delimiter = ","), "gzip|read|compressed")
    expect_false(dir.exists(failed))
})
