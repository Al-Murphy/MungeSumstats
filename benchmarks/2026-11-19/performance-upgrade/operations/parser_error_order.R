# Additional installed public probe. Run each role in a fresh process and directory.
args <- commandArgs(TRUE)
stopifnot(length(args) == 5L, args[[5L]] %in% c("width", "precedence"))
.libPaths(c(args[[1L]], args[[2L]]), include.site = FALSE)
library(MungeSumstats, lib.loc = args[[1L]])
root <- args[[4L]]
stopifnot(!file.exists(root), dir.create(root))

hdr <- "CHR\tBP\n"
ragged <- "S\t2\tx\n"
valid <- "S\t2\n"
late_quote <- 'S\tq"\n'
width <- "Input row has an unexpected number of fields near source row "
bad_quote <- "Malformed quote in unquoted field"
simple <- c("simpleError", "error", "condition")
native <- c("Rcpp::exception", "C++Error", "error", "condition")
records <- list()

for (format in 2:3) {
    handle <- open_reference(file.path(args[[3L]], paste0("scientific_", format)))
    probe <- function(name, bytes, rows, budget, want, classes, cap = 128L) {
        id <- paste0(name, "-f", format)
        input <- file.path(root, paste0(id, ".tsv"))
        dest <- file.path(root, paste0(id, "-edition"))
        writeBin(bytes, input)
        public_call <- quote(harmonise_file(input, handle, dest, batch_size = rows,
                                            batch_bytes = budget, max_record_bytes = cap))
        err <- tryCatch({ eval(public_call); NULL }, error = identity)
        msg <- if (is.null(err)) "" else conditionMessage(err)
        published <- file.exists(dest)
        # Width belongs to harmonise_file's caller. Lexical errors come from the native decoder.
        width_call <- if (startsWith(msg, width)) identical(conditionCall(err), public_call) else NA
        ok <- identical(class(err), classes) && identical(msg, want) && !published &&
              (!startsWith(want, width) || isTRUE(width_call))
        records[[id]] <<- data.frame(case = id, class = paste(class(err), collapse = ","),
            error = msg, width_call = width_call, published = published, pass = ok)
    }
    if (args[[5L]] == "width") {
        probe("first", charToRaw(paste0(hdr, ragged)), 2L, 64L, paste0(width, 1L), simple)
        # Row four fails in the batch beginning at row three.
        probe("ordinal", charToRaw(paste0(hdr, valid, valid, valid, ragged)), 2L, 64L,
              paste0(width, 3L), simple)
    } else {
        mixed <- charToRaw(paste0(hdr, ragged, late_quote))
        probe("same-block", mixed, 1L, 64L, bad_quote, native)
        # Thirteen bytes end immediately after the ragged row. A full batch
        # stops there; a partially filled batch feeds the next block first.
        probe("across-full", mixed, 1L, 13L, paste0(width, 1L), simple)
        probe("across-partial", mixed, 2L, 13L, bad_quote, native)
        probe("later-nul", c(charToRaw(paste0(hdr, ragged)), as.raw(0L)),
              1L, 64L, "NUL byte in tabular input", native)
        probe("later-limit", charToRaw(paste0(hdr, ragged, "12345678")), 1L, 64L,
              "Input record exceeds max_record_bytes", native, cap = 7L)
        eof_quote <- charToRaw(paste0(hdr, ragged, 'S\t"q'))
        probe("eof-full", eof_quote, 1L, 64L, paste0(width, 1L), simple)
        probe("eof-partial", eof_quote, 2L, 64L, "Unterminated quoted field", native)
        probe("duplicate-header", charToRaw(paste0("BP\tBP\n", late_quote)), 1L, 64L,
              bad_quote, native)
    }
}
results <- do.call(rbind, records)
utils::write.table(results, file.path(root, "errors.tsv"), sep = "\t", row.names = FALSE,
                   quote = TRUE, qmethod = "double")
stopifnot(all(results$pass))
