#' @useDynLib MungeSumstats, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Prepare a portable, preservation-oriented reference
#'
#' Builds an immutable disk index without loading all variants or FASTA bases in
#' memory. This new interface does not change [format_sumstats()].
#' @param fasta Path to a plain, fixed-line-width FASTA. Compressed FASTA is not
#'   supported. Contig names must agree exactly with the reference and input.
#' @param variants Path to a TSV with exactly these columns in order:
#'   `chr,pos,ref,alt,id,aliases`. Positions are 1-based; ALT may contain comma
#'   separated alleles and aliases are semicolon separated. `.` denotes missing.
#' @param path New, version-specific output directory. Existing paths are refused.
#' @param build Explicit coordinate build label, e.g. `GRCh38`.
#' @param version Explicit reference release/version label.
#' @param batch_size Maximum index sort buffer rows. Merge fan-in is 32.
#' @return An `mss_reference` containing build, version, paths and MD5 checksums.
#' @details Each valid REF/ALT pair is independently trimmed and left aligned
#'   against FASTA. Every rsID/alias can resolve to multiple candidates. Invalid or
#'   unsupported reference alleles appear in `rejected.tsv`; the exact reference
#'   TSV is retained as `source.tsv`. FASTA remains at its declared absolute path.
#'   Metadata and checksums identify that external dependency; do not mutate it.
#'   Memory includes at most the sort buffer, 32 merge rows, contig metadata and
#'   individual source lines, plus logarithmically many temporary run paths. A
#'   single unusually large line can exceed the nominal row-buffer bound.
#' @export
prepare_reference <- function(fasta, variants, path, build, version,
                              batch_size = 100000L) {
    if (missing(build) || length(build) != 1L || is.na(build) || !nzchar(build))
        stop("An explicit build is required")
    if (missing(version) || length(version) != 1L || is.na(version) || !nzchar(version))
        stop("An explicit version is required")
    batch_size <- .mss_positive_integer(batch_size, "batch_size")
    fasta <- normalizePath(fasta, mustWork = TRUE)
    variants <- normalizePath(variants, mustWork = TRUE)
    if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("Invalid path")
    if (file.exists(path)) stop("Reference path already exists; choose a new version path")
    parent <- normalizePath(dirname(path), mustWork = TRUE)
    path <- file.path(parent, basename(path))
    staging <- tempfile(paste0(".", basename(path), "-pending-"), tmpdir = parent)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    input_hashes <- unname(tools::md5sum(c(fasta, variants)))
    mss_prepare_cpp(fasta, variants, staging, batch_size)
    if (!file.copy(variants, file.path(staging, "source.tsv"))) stop("Cannot preserve reference TSV")
    if (!identical(input_hashes, unname(tools::md5sum(c(fasta, variants)))) ||
        !identical(input_hashes[[2L]], unname(tools::md5sum(file.path(staging, "source.tsv")))))
        stop("Reference input changed during indexing; no reference was published")
    filenames <- c("fasta.fai", "keys", "offsets", "records", "rejected.tsv", "source.tsv", "core-version")
    files <- c(fasta = fasta, stats::setNames(file.path(staging, filenames), filenames))
    checksums <- tools::md5sum(files)
    names(checksums) <- names(files)
    if (!identical(checksums[["fasta"]], input_hashes[[1L]]) ||
        !identical(checksums[["source.tsv"]], input_hashes[[2L]]))
        stop("Reference input changed while hashing; no reference was published")
    sizes <- file.info(files)$size
    names(sizes) <- names(files)
    manifest <- list(schema_version = 1L, core_version = "mss-preservation-2",
                     build = build, version = version, fasta = fasta,
                     checksum_algorithm = "MD5", checksums = checksums,
                     sizes = sizes, files = filenames,
                     created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     input_variants = variants,
                     projection = "not_attempted")
    saveRDS(manifest, file.path(staging, "manifest.rds"), version = 2)
    if (file.exists(path) || !file.rename(staging, path)) stop("Cannot publish immutable reference path")
    .mss_read_reference(path)
}

.mss_positive_integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
        x < 1 || x > .Machine$integer.max || x != floor(x))
        stop(name, " must be a positive integer")
    as.integer(x)
}

.mss_read_reference <- function(reference) {
    if (inherits(reference, "mss_reference_handle")) return(reference)
    path <- if (inherits(reference, "mss_reference")) reference$path else reference
    if (!is.character(path) || length(path) != 1L) stop("reference must be an mss_reference or its path")
    path <- normalizePath(path, mustWork = TRUE)
    manifest <- readRDS(file.path(path, "manifest.rds"))
    if (!identical(manifest$schema_version, 1L) ||
        !identical(manifest$core_version, "mss-preservation-2")) stop("Unsupported reference schema")
    manifest$path <- path
    class(manifest) <- "mss_reference"
    manifest
}

.mss_verify_reference <- function(reference, verify) {
    files <- c(fasta = reference$fasta,
               stats::setNames(file.path(reference$path, reference$files), reference$files))
    info <- file.info(files)
    if (inherits(reference, "mss_reference_handle")) {
        if (!identical(info[, c("size", "mtime")], reference$.file_signature))
            stop("Reference checksum/stat mismatch; reopen a restored immutable reference")
        return(invisible(TRUE))
    }
    sizes <- info$size
    if (anyNA(sizes) || !identical(as.numeric(sizes), as.numeric(reference$sizes)))
        stop("Reference checksum/size mismatch; restore the immutable reference")
    if (verify) {
        checksums <- tools::md5sum(files)
        if (!identical(unname(checksums), unname(reference$checksums)))
            stop("Reference checksum mismatch; restore the immutable reference")
    }
    invisible(TRUE)
}

#' Open a verified reusable disk reference handle
#'
#' @param reference An `mss_reference` or its directory path.
#' @return A live `mss_reference_handle` accepted by both harmonisation APIs.
#' @details MD5 verification occurs once when opening; subsequent uses verify
#'   sizes and modification times. The handle owns native disk file streams and
#'   contig metadata, not a reference-variant or FASTA sequence cache. Keep the
#'   reference files immutable while the handle is live. Handles are process-local
#'   and must be reopened after serialization or in another R process.
#' @export
open_reference <- function(reference) {
    reference <- .mss_read_reference(reference)
    .mss_verify_reference(reference, TRUE)
    if (inherits(reference, "mss_reference_handle")) return(reference)
    files <- c(fasta = reference$fasta,
               stats::setNames(file.path(reference$path, reference$files), reference$files))
    signature <- file.info(files)[, c("size", "mtime")]
    reference$.native <- mss_open_reference_cpp(reference$path, reference$fasta)
    if (!identical(signature, file.info(files)[, c("size", "mtime")]))
        stop("Reference changed while opening its native handle")
    reference$.file_signature <- signature
    class(reference) <- c("mss_reference_handle", "mss_reference")
    reference
}

#' Harmonise summary statistics while retaining every source row
#'
#' This separate API preserves the input exactly and returns a row-aligned,
#' nullable normalized table. It never calls the legacy filtering pipeline.
#' @param sumstats A data.frame, data.table, or previous `mss_harmonised` result.
#' @param reference An `mss_reference` or reference directory path.
#' @param columns Named list mapping `chr,pos,id,other,effect_allele,effect,eaf`
#'   to input column names. Defaults use `CHR,BP,SNP,A1,A2,BETA,FRQ`; `A2` is the
#'   effect allele. Missing columns are allowed and remain unknown.
#' @param build Explicit input build. Defaults to the declared reference build;
#'   a different build is rejected. No liftover is performed.
#' @param frequency_type `EAF`, `MAF`, or `unknown`. Only an explicit EAF can be
#'   transformed. A swapped EAF is complemented only at a known biallelic site.
#' @param effect_type `beta`, `log_odds`, or `none`. Only signed additive effects
#'   can be negated; odds ratios must first be explicitly converted by the caller.
#' @param batch_size Maximum input rows sent to each native call.
#' @param candidate_limit Maximum raw index candidates per disk lookup; overflow stops
#'   processing rather than silently selecting a candidate.
#' @param verify Whether to verify all reference MD5 hashes (default `TRUE`).
#'   `FALSE` checks sizes only and is appropriate only for already trusted,
#'   immutable files; full verification reads the reference files sequentially.
#'   A live `open_reference()` handle reuses its initial verification and checks
#'   file sizes and modification times on subsequent calls.
#' @return An `mss_harmonised` list: exact `source`, row-aligned `normalized`, and
#'   `provenance`. No source measurements, identifiers, optional fields or raw
#'   tokens are overwritten. Unresolved normalized measurements are `NA`.
#' @details Recognized outcomes include `match`, `swap`, `complement`,
#'   `swapped_complement`, `alt_alt_contrast`, `ambiguous`, `unmapped`,
#'   `allele_mismatch`, `unsupported_allele`, `invalid_position`, and
#'   `coordinate_id_conflict`. Palindromic contrasts that admit both effect
#'   directions remain ambiguous. ALT-ALT contrasts are retained without
#'   inventing a reference-anchored effect. All GRCh38 projection fields remain
#'   missing even when native coordinates are GRCh38: projection was not attempted.
#'   Symbolic alleles and unknown IDs are preserved in `source`. Input and output
#'   tables remain resident in R; batching bounds native work, not total R memory.
#' @export
harmonise_sumstats <- function(sumstats, reference,
                              columns = list(chr = "CHR", pos = "BP", id = "SNP",
                                             other = "A1", effect_allele = "A2",
                                             effect = "BETA", eaf = "FRQ"),
                              build = NULL, frequency_type = c("unknown", "EAF", "MAF"),
                              effect_type = c("beta", "log_odds", "none"),
                              batch_size = 100000L, candidate_limit = 10000L,
                              verify = TRUE) {
    if (inherits(sumstats, "mss_harmonised")) sumstats <- sumstats$source
    if (!is.data.frame(sumstats)) stop("sumstats must be a data.frame or data.table")
    frequency_type <- match.arg(frequency_type)
    effect_type <- match.arg(effect_type)
    batch_size <- .mss_positive_integer(batch_size, "batch_size")
    candidate_limit <- .mss_positive_integer(candidate_limit, "candidate_limit")
    if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) stop("verify must be TRUE or FALSE")
    reference <- .mss_read_reference(reference)
    if (is.null(build)) build <- reference$build
    if (!identical(build, reference$build)) stop("Input build differs from reference build; projection is not implemented")
    .mss_verify_reference(reference, verify)
    required <- c("chr", "pos", "id", "other", "effect_allele", "effect", "eaf")
    if (!is.list(columns) || length(names(columns)) != length(columns) ||
        any(!nzchar(names(columns))) || any(!names(columns) %in% required) || anyDuplicated(names(columns)))
        stop("columns must be a named list using chr,pos,id,other,effect_allele,effect,eaf")
    n <- nrow(sumstats)
    extract <- function(field, rows) {
        name <- columns[[field]]
        if (is.null(name)) return(rep(NA_character_, length(rows)))
        if (!is.character(name) || length(name) != 1L || is.na(name)) stop("Column mappings must be single column names or NULL")
        if (!name %in% names(sumstats)) return(rep(NA_character_, length(rows)))
        if (anyDuplicated(names(sumstats)) && sum(names(sumstats) == name) > 1L) stop("Mapped column name is duplicated: ", name)
        as.character(sumstats[[name]][rows])
    }
    parse_number <- function(value) suppressWarnings(as.numeric(value))
    batches <- vector("list", if (n) ceiling(n / batch_size) else 1L)
    starts <- if (n) seq.int(1, n, by = batch_size) else 1L
    for (i in seq_along(starts)) {
        rows <- if (n) seq.int(starts[[i]], min(n, starts[[i]] + batch_size - 1)) else integer()
        raw_pos <- extract("pos", rows)
        pos <- parse_number(raw_pos)
        # Distinguish an unknown position from a malformed token.
        invalid_pos <- !is.na(raw_pos) & !raw_pos %in% c("", ".", "NA") & is.na(pos)
        pos[invalid_pos] <- -1
        effect <- if (effect_type == "none") rep(NA_real_, length(rows)) else parse_number(extract("effect", rows))
        batch <- mss_harmonise_cpp(reference$path, reference$fasta,
                                  extract("chr", rows), pos, extract("id", rows),
                                  extract("other", rows), extract("effect_allele", rows),
                                  effect, parse_number(extract("eaf", rows)),
                                  frequency_type, candidate_limit,
                                  if (inherits(reference, "mss_reference_handle")) reference$.native else NULL)
        batch$normalized_build <- rep(build, length(rows))
        batch$source_row <- as.double(rows)
        batch$grch38_chr <- rep(NA_character_, length(rows))
        batch$grch38_pos <- rep(NA_real_, length(rows))
        batch$grch38_ref <- rep(NA_character_, length(rows))
        batch$grch38_alt <- rep(NA_character_, length(rows))
        batches[[i]] <- batch
    }
    normalized <- do.call(rbind, batches)
    rownames(normalized) <- NULL
    provenance <- list(core_version = reference$core_version,
                       reference_version = reference$version, build = build,
                       reference_path = reference$path,
                       reference_checksums = reference$checksums,
                       frequency_type = frequency_type, effect_type = effect_type,
                       columns = columns, reference_verified = verify,
                       projection = "not_attempted", schema_version = 1L)
    structure(list(source = sumstats, normalized = normalized, provenance = provenance),
              class = "mss_harmonised")
}

#' Stream a preservation-oriented harmonisation edition to disk
#'
#' @param input Plain or gzip-compressed tabular input file.
#' @param reference An `mss_reference` or reference directory.
#' @param path New immutable output directory.
#' @param columns,build,frequency_type,effect_type,candidate_limit See
#'   [harmonise_sumstats()].
#' @param delimiter A single tab (default) or comma. Input is explicitly
#'   unquoted, one record per line; quoted/multiline fields are rejected.
#' @param batch_size Maximum source lines resident in each batch.
#' @return A manifest with output path, source row count and status counts.
#' @details The original input bytes are copied to `source-artifact`; decompressed
#'   raw tokens are retained in `source.tsv` (or `source.csv`). `normalized.tsv`
#'   carries a 1-based `source_row` association and nullable native fields.
#'   Scientific numeric output uses 17 significant decimal digits. `audit.tsv`
#'   counts every status. `manifest.rds` records checksums and reference identity.
#'   Portable quoted TSV receipts are also written: `metadata.tsv` (key/value),
#'   `columns.tsv` (field/source_column), `reference_checksums.tsv` (name/md5),
#'   and `checksums.tsv` (filename/md5). The checksum list excludes itself and
#'   the RDS manifest. `normalized_build` explicitly labels native coordinates.
#'   No complete input/output table is loaded in R. Memory scales with batch size,
#'   individual line lengths and the native candidate cap; source lines are not
#'   artificially truncated. Reference hashes are verified once before streaming.
#'   A directory is published only after every batch and final audit succeeds.
#'   Malformed field counts, quoted fields and read warnings (including missing
#'   final newline or detected truncated gzip input) fail without publication.
#'   Syntactically valid but incomplete files cannot be distinguished from intended
#'   smaller inputs without an externally supplied acquisition receipt.
#' @export
harmonise_file <- function(input, reference, path,
                           columns = list(chr = "CHR", pos = "BP", id = "SNP",
                                          other = "A1", effect_allele = "A2",
                                          effect = "BETA", eaf = "FRQ"),
                           build = NULL, frequency_type = c("unknown", "EAF", "MAF"),
                           effect_type = c("beta", "log_odds", "none"),
                           delimiter = "\t", batch_size = 100000L,
                           candidate_limit = 10000L) {
    batch_size <- .mss_positive_integer(batch_size, "batch_size")
    candidate_limit <- .mss_positive_integer(candidate_limit, "candidate_limit")
    frequency_type <- match.arg(frequency_type)
    effect_type <- match.arg(effect_type)
    if (!identical(delimiter, "\t") && !identical(delimiter, ",")) stop("delimiter must be tab or comma")
    input <- normalizePath(input, mustWork = TRUE)
    reference <- .mss_read_reference(reference)
    if (is.null(build)) build <- reference$build
    if (!identical(build, reference$build)) stop("Input build differs from reference build")
    reference <- open_reference(reference)
    if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("Invalid path")
    if (file.exists(path)) stop("Output path already exists; choose a new edition path")
    parent <- normalizePath(dirname(path), mustWork = TRUE)
    path <- file.path(parent, basename(path))
    staging <- tempfile(paste0(".", basename(path), "-pending-"), tmpdir = parent)
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    before <- unname(tools::md5sum(input))
    snapshot <- file.path(staging, "source-artifact")
    if (!file.copy(input, snapshot) ||
        !identical(before, unname(tools::md5sum(snapshot))) ||
        !identical(before, unname(tools::md5sum(input)))) stop("Input changed during snapshot")
    reader <- gzfile(snapshot, open = "rt")
    on.exit(close(reader), add = TRUE)
    read_batch <- function(n) withCallingHandlers(readLines(reader, n = n, warn = TRUE),
                                                 warning = function(w) stop("Input read/line failure: ", conditionMessage(w)))
    header <- read_batch(1L)
    if (length(header) != 1L || !nzchar(header)) stop("Missing input header")
    split_line <- function(lines) {
        # Appending a sentinel preserves empty trailing fields with strsplit.
        lapply(strsplit(paste0(lines, delimiter, "__MSS_END__"), delimiter, fixed = TRUE),
               function(fields) fields[-length(fields)])
    }
    if (grepl('"', header, fixed = TRUE)) stop("Quoted fields are unsupported")
    header_names <- split_line(header)[[1L]]
    if (any(!nzchar(header_names)) || anyDuplicated(header_names)) stop("Header names must be nonempty and unique")
    source_name <- if (delimiter == "\t") "source.tsv" else "source.csv"
    source_writer <- file(file.path(staging, source_name), "wt")
    on.exit(close(source_writer), add = TRUE)
    normalized_writer <- file(file.path(staging, "normalized.tsv"), "wt")
    on.exit(close(normalized_writer), add = TRUE)
    writeLines(header, source_writer)
    total <- 0
    counts <- numeric()
    first <- TRUE
    write_normalized <- function(table) {
        # write.table's numeric conversion can round to 15 digits; serialize
        # float64 measurements explicitly before writing delimited text.
        encoded <- lapply(table, function(value) {
            if (is.numeric(value)) {
                result <- sprintf("%.17g", value)
                result[is.na(value)] <- NA_character_
                result
            } else as.character(value)
        })
        encoded <- as.data.frame(encoded, stringsAsFactors = FALSE, check.names = FALSE)
        utils::write.table(encoded, normalized_writer, sep = "\t", quote = TRUE,
                           row.names = FALSE, col.names = first, na = "NA", qmethod = "double")
        first <<- FALSE
    }
    repeat {
        lines <- read_batch(batch_size)
        if (!length(lines)) break
        if (any(grepl('"', lines, fixed = TRUE))) stop("Quoted fields are unsupported")
        fields <- split_line(lines)
        if (any(lengths(fields) != length(header_names))) stop("Input row has an unexpected number of fields near source row ", total + 1)
        source_batch <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE, check.names = FALSE)
        names(source_batch) <- header_names
        result <- harmonise_sumstats(source_batch, reference, columns = columns, build = build,
                                     frequency_type = frequency_type, effect_type = effect_type,
                                     batch_size = batch_size, candidate_limit = candidate_limit,
                                     verify = FALSE)$normalized
        result$source_row <- result$source_row + total
        write_normalized(result)
        writeLines(lines, source_writer)
        batch_counts <- table(result$status)
        for (status in names(batch_counts)) {
            if (!status %in% names(counts)) counts[[status]] <- 0
            counts[[status]] <- counts[[status]] + as.double(batch_counts[[status]])
        }
        total <- total + length(lines)
    }
    if (first) {
        empty <- as.data.frame(stats::setNames(rep(list(character()), length(header_names)), header_names),
                               stringsAsFactors = FALSE, check.names = FALSE)
        write_normalized(harmonise_sumstats(empty, reference, columns = columns, build = build,
                                           frequency_type = frequency_type, effect_type = effect_type,
                                           candidate_limit = candidate_limit, verify = FALSE)$normalized)
    }
    flush(source_writer); flush(normalized_writer)
    if (sum(counts) != total) stop("Source/native audit row count mismatch")
    audit <- data.frame(status = names(counts), count = as.double(counts))
    utils::write.table(audit, file.path(staging, "audit.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE)
    metadata <- c(schema_version = "1", core_version = reference$core_version,
                  reference_build = build, reference_version = reference$version,
                  reference_path = reference$path, source_rows = sprintf("%.17g", total),
                  normalized_rows = sprintf("%.17g", total), input_checksum = before,
                  checksum_algorithm = "MD5", frequency_type = frequency_type,
                  effect_type = effect_type, projection = "not_attempted",
                  delimiter = if (delimiter == "\t") "tab" else "comma",
                  source_file = source_name)
    utils::write.table(data.frame(key = names(metadata), value = unname(metadata)),
                       file.path(staging, "metadata.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE, qmethod = "double")
    column_names <- vapply(columns, function(column) if (is.null(column)) NA_character_ else column,
                           character(1))
    utils::write.table(data.frame(field = names(columns), source_column = column_names),
                       file.path(staging, "columns.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE, qmethod = "double")
    utils::write.table(data.frame(name = names(reference$checksums), md5 = unname(reference$checksums)),
                       file.path(staging, "reference_checksums.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE, qmethod = "double")
    filenames <- c("source-artifact", source_name, "normalized.tsv", "audit.tsv",
                   "metadata.tsv", "columns.tsv", "reference_checksums.tsv")
    checksums <- tools::md5sum(file.path(staging, filenames))
    names(checksums) <- filenames
    utils::write.table(data.frame(filename = filenames, md5 = unname(checksums)),
                       file.path(staging, "checksums.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE, qmethod = "double")
    manifest <- list(schema_version = 1L, core_version = reference$core_version,
                     reference_build = build, reference_version = reference$version,
                     reference_checksums = reference$checksums, reference_path = reference$path,
                     source_rows = total, normalized_rows = total,
                     status_counts = counts, input_checksum = before,
                     checksums = checksums, checksum_algorithm = "MD5",
                     columns = columns, frequency_type = frequency_type,
                     effect_type = effect_type, projection = "not_attempted",
                     delimiter = delimiter,
                     created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    saveRDS(manifest, file.path(staging, "manifest.rds"), version = 2)
    # Close output connections before rename, which is required on Windows.
    close(normalized_writer);close(source_writer);close(reader)
    on.exit(NULL, add = FALSE)
    if (file.exists(path) || !file.rename(staging, path)) {
        unlink(staging, recursive = TRUE)
        stop("Cannot publish immutable output path")
    }
    manifest$path <- path
    manifest
}
