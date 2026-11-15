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
#' @param contig_aliases Optional named character vector from authoritative
#'   source labels to exact FASTA contig names. Many aliases may name one contig.
#'   Existing FASTA names cannot be redirected to another contig. Empty mapping
#'   retains strict-name behavior; no aliases or build changes are inferred.
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
                              batch_size = 100000L, contig_aliases = character()) {
    if (missing(build) || length(build) != 1L || is.na(build) || !nzchar(build))
        stop("An explicit build is required")
    if (missing(version) || length(version) != 1L || is.na(version) || !nzchar(version))
        stop("An explicit version is required")
    batch_size <- .mss_positive_integer(batch_size, "batch_size")
    contig_aliases <- .mss_validate_contig_aliases(contig_aliases)
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
    if (length(contig_aliases)) {
        fai <- utils::read.delim(file.path(staging, "fasta.fai"), header = FALSE,
                                colClasses = "character", comment.char = "")
        .mss_validate_contig_aliases(contig_aliases, fai[[1L]])
        utils::write.table(data.frame(alias = names(contig_aliases),
                                      canonical_contig = unname(contig_aliases)),
                           file.path(staging, "contig_aliases.tsv"), sep = "\t",
                           quote = TRUE, row.names = FALSE, qmethod = "double")
    }
    if (!file.copy(variants, file.path(staging, "source.tsv"))) stop("Cannot preserve reference TSV")
    if (!identical(input_hashes, unname(tools::md5sum(c(fasta, variants)))) ||
        !identical(input_hashes[[2L]], unname(tools::md5sum(file.path(staging, "source.tsv")))))
        stop("Reference input changed during indexing; no reference was published")
    filenames <- c("fasta.fai", "keys", "offsets", "records", "rejected.tsv", "source.tsv", "core-version")
    if (length(contig_aliases)) filenames <- c(filenames, "contig_aliases.tsv")
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

.mss_validate_contig_aliases <- function(x, contigs = NULL) {
    if (!is.character(x) || anyNA(x) || any(!nzchar(x)) ||
        (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
                       any(!nzchar(names(x))) || anyDuplicated(names(x)))))
        stop("Contig aliases must be a uniquely named character vector")
    if (!length(x)) return(character())
    if (length(x) && any(toupper(names(x)) %in% c("NA", "N/A", "NAN", "NULL", ".", "-")))
        stop("Contig aliases cannot map missing-value labels")
    if (!is.null(contigs)) {
        if (any(!x %in% contigs)) stop("Contig alias target is absent from FASTA")
        canonical <- names(x) %in% contigs
        if (any(names(x)[canonical] != x[canonical]))
            stop("Contig aliases cannot redirect an existing FASTA name")
    }
    x[order(names(x))]
}

.mss_alias_contigs <- function(value, reference) {
    aliases <- reference$contig_aliases
    if (!length(aliases)) return(value)
    matched <- match(value, names(aliases))
    known <- !is.na(matched)
    value[known] <- unname(aliases[matched[known]])
    value
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
    manifest$contig_aliases <- character()
    if ("contig_aliases.tsv" %in% manifest$files) {
        aliases <- utils::read.delim(file.path(path, "contig_aliases.tsv"),
                                    colClasses = "character", comment.char = "",
                                    check.names = FALSE)
        if (!identical(names(aliases), c("alias", "canonical_contig")))
            stop("Invalid contig alias table")
        manifest$contig_aliases <- .mss_validate_contig_aliases(
            stats::setNames(aliases$canonical_contig, aliases$alias))
    }
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
#'   contig metadata, plus at most 28 MiB of reference/FASTA byte pages (and
#'   bounded page metadata), independently of total reference size. Keep the
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
#'   to input column names. Optional `entity_type,entity_id` preserve specialist
#'   statistical entities. Types are `small_variant`, `symbolic_sv`, `breakend`,
#'   `hla`, `repeat`, `gene_burden`; only SV/BND syntax is automatically typed.
#'   Specialist rows retain identity and source tokens with status
#'   `specialist_preserved`, without small-variant normalized effects. Defaults use `CHR,BP,SNP,A1,A2,BETA,FRQ`; `A2` is the
#'   effect allele. Missing columns are allowed and remain unknown. Optional
#'   `source_ref,source_alt` mappings declare trusted source reference/alternate
#'   alleles, together with either a `source_strand` column or the static argument.
#' @param source_strand Optional explicit file/table-level `"+"`, `"-"` or
#'   `"unknown"` declaration for trusted source REF/ALT. Cannot be combined with
#'   a `columns$source_strand` mapping. Coordinates always name the forward
#'   reference's leftmost 1-based position. REF/ALT and measured allele strings
#'   must share the declared strand. No declaration is inferred from A1/A2.
#' @param build Explicit input build. Defaults to the declared reference build;
#'   a different build is rejected. Use `projection` for an explicit target build.
#' @param frequency_type `EAF`, `MAF`, or `unknown`. Only an explicit EAF can be
#'   transformed. A swapped EAF is complemented only at a known biallelic site.
#' @param effect_type `beta`, `log_odds`, or `none`. Only signed additive effects
#'   can be negated; odds ratios must first be explicitly converted by the caller.
#' @param batch_size Maximum input rows sent to each native call.
#' @param candidate_limit Maximum raw index candidates per disk lookup; overflow stops
#'   processing rather than silently selecting a candidate.
#' @param projection Optional handle from [prepare_projection()].
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
#'   missing unless an explicit `projection` is supplied. Projected coordinates,
#'   alleles and effects have separate columns; native fields are never replaced.
#'   Symbolic alleles and unknown IDs are preserved in `source`. Input and output
#'   tables remain resident in R; batching bounds native work, not total R memory.
#'   With an explicit source REF/ALT/strand contract, source REF must validate
#'   against FASTA before normalization. Unknown strand, allele/ID conflict or
#'   source REF mismatch remains unresolved. Valid sequence absent from the
#'   catalog is labeled `sequence_only`, without inventing an identifier.
#'   Generic indels whose opposite event directions both validate against FASTA
#'   require a resolving known identifier; otherwise they remain ambiguous.
#' @export
harmonise_sumstats <- function(sumstats, reference,
                              columns = list(chr = "CHR", pos = "BP", id = "SNP",
                                             other = "A1", effect_allele = "A2",
                                             effect = "BETA", eaf = "FRQ"),
                              build = NULL, frequency_type = c("unknown", "EAF", "MAF"),
                              effect_type = c("beta", "log_odds", "none"),
                              batch_size = 100000L, candidate_limit = 10000L,
                              verify = TRUE, source_strand = NULL, projection = NULL) {
    if (inherits(sumstats, "mss_harmonised")) sumstats <- sumstats$source
    if (!is.data.frame(sumstats)) stop("sumstats must be a data.frame or data.table")
    frequency_type <- match.arg(frequency_type)
    effect_type <- match.arg(effect_type)
    batch_size <- .mss_positive_integer(batch_size, "batch_size")
    candidate_limit <- .mss_positive_integer(candidate_limit, "candidate_limit")
    if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) stop("verify must be TRUE or FALSE")
    reference <- .mss_read_reference(reference)
    if (is.null(build)) build <- reference$build
    if (!identical(build, reference$build)) stop("Input build differs from native reference build; supply projection separately")
    .mss_verify_reference(reference, verify)
    required <- c("chr", "pos", "id", "other", "effect_allele", "effect", "eaf",
                  "source_ref", "source_alt", "source_strand", "entity_type", "entity_id")
    if (!is.list(columns) || length(names(columns)) != length(columns) ||
        any(!nzchar(names(columns))) || any(!names(columns) %in% required) || anyDuplicated(names(columns)))
        stop("columns must be a named list using chr,pos,id,other,effect_allele,effect,eaf")
    source_contract <- any(c("source_ref", "source_alt", "source_strand") %in% names(columns)) || !is.null(source_strand)
    if (source_contract) {
        if (!all(c("source_ref", "source_alt") %in% names(columns)) ||
            (is.null(source_strand) && !"source_strand" %in% names(columns)))
            stop("source_ref, source_alt and a source strand declaration must be supplied together")
        if (!is.null(source_strand) && "source_strand" %in% names(columns))
            stop("Choose either a source_strand column or a file-level source_strand declaration")
        if (!is.null(source_strand) && (!is.character(source_strand) || length(source_strand) != 1L ||
                                      is.na(source_strand) || !source_strand %in% c("+", "-", "unknown")))
            stop("source_strand must be '+', '-' or 'unknown'")
    }
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
        # Explicit statistical entities never acquire a small-variant effect.
        entity <- tolower(extract("entity_type", rows))
        entity[is.na(entity) | entity %in% c("", ".", "na")] <- "small_variant"
        alleles <- paste(extract("other", rows), extract("effect_allele", rows))
        entity[entity == "small_variant" & grepl("<[^>]+>", alleles)] <- "symbolic_sv"
        entity[entity == "small_variant" & grepl("[\\[\\]]", alleles, perl = TRUE)] <- "breakend"
        supported <- c("small_variant", "symbolic_sv", "breakend", "hla", "repeat", "gene_burden")
        if (any(!entity %in% supported)) stop("Unknown entity_type; use small_variant, symbolic_sv, breakend, hla, repeat or gene_burden")
        specialist <- entity != "small_variant"
        native_extract <- function(field, rows) {
            value <- extract(field, rows)
            value[specialist] <- NA_character_
            value
        }
        raw_pos <- native_extract("pos", rows)
        pos <- parse_number(raw_pos)
        # Distinguish an unknown position from a malformed token.
        invalid_pos <- !is.na(raw_pos) & !raw_pos %in% c("", ".", "NA") & is.na(pos)
        pos[invalid_pos] <- -1
        effect <- if (effect_type == "none") rep(NA_real_, length(rows)) else parse_number(native_extract("effect", rows))
        batch <- mss_harmonise_cpp(reference$path, reference$fasta,
                                  .mss_alias_contigs(native_extract("chr", rows), reference), pos, native_extract("id", rows),
                                  native_extract("other", rows), native_extract("effect_allele", rows),
                                  effect, parse_number(native_extract("eaf", rows)),
                                  frequency_type, candidate_limit,
                                  if (inherits(reference, "mss_reference_handle")) reference$.native else NULL,
                                  if (source_contract) native_extract("source_ref", rows) else NULL,
                                  if (source_contract) native_extract("source_alt", rows) else NULL,
                                  if (source_contract) {
                                      if (is.null(source_strand)) native_extract("source_strand", rows) else rep(source_strand, length(rows))
                                  } else NULL)
        batch$status[specialist] <- "specialist_preserved"
        for (field in setdiff(names(batch), c("status", "candidate_count"))) batch[[field]][specialist] <- NA
        batch$candidate_count[specialist] <- 0
        batch$entity_type <- entity
        batch$entity_id <- extract("entity_id", rows)
        absent_id <- is.na(batch$entity_id)
        batch$entity_id[absent_id] <- extract("id", rows)[absent_id]
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
    if (!is.null(projection)) normalized <- .mss_project(normalized, reference, projection, candidate_limit)
    provenance <- list(core_version = reference$core_version,
                       implementation_version = "mss-preservation-4",
                       reference_version = reference$version, build = build,
                       reference_path = reference$path,
                       reference_checksums = reference$checksums,
                       contig_mapping_implementation = if (length(reference$contig_aliases)) "mss-contigs-1" else "not_applied",
                       frequency_type = frequency_type, effect_type = effect_type,
                       columns = columns, reference_verified = verify,
                       projection = if (is.null(projection)) "not_attempted" else projection$receipt, schema_version = 1L,
                       source_strand = source_strand)
    structure(list(source = sumstats, normalized = normalized, provenance = provenance),
              class = "mss_harmonised")
}

#' Stream a preservation-oriented harmonisation edition to disk
#'
#' @param input Plain or gzip-compressed tabular input file.
#' @param reference An `mss_reference` or reference directory.
#' @param path New immutable output directory.
#' @param columns,build,frequency_type,effect_type,candidate_limit,source_strand See
#'   [harmonise_sumstats()].
#' @param delimiter A single tab (default) or comma. Double-quoted fields, doubled
#'   escaped quotes, embedded delimiters/newlines, CRLF and final records without
#'   newlines are supported. Quotes in unquoted fields are syntax errors.
#' @param batch_size Maximum source records resident in each batch.
#' @param batch_bytes Adaptive decoded-token byte budget (default 8 MiB). A
#'   single larger record is allowed up to `max_record_bytes`.
#' @param max_record_bytes Maximum encoded bytes per record (default 64 MiB).
#'   Oversized records fail without publication; increase explicitly if needed.
#' @param acquisition_receipt Optional named list of expected `rows`, original
#'   compressed or plain `bytes`, and/or original file `md5`. Mismatch aborts.
#' @param projection Optional handle from [prepare_projection()].
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
#'   byte/record limits and the native candidate cap. No record is truncated. Reference hashes are verified once before streaming.
#'   A directory is published only after every batch and final audit succeeds.
#'   Malformed quotes/field counts and read warnings (including detected truncated
#'   gzip input) fail without publication. A missing final LF is valid syntax.
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
                           candidate_limit = 10000L, source_strand = NULL,
                           projection = NULL, batch_bytes = 8L * 1024L * 1024L,
                           max_record_bytes = 64L * 1024L * 1024L,
                           acquisition_receipt = NULL) {
    batch_size <- .mss_positive_integer(batch_size, "batch_size")
    candidate_limit <- .mss_positive_integer(candidate_limit, "candidate_limit")
    frequency_type <- match.arg(frequency_type)
    effect_type <- match.arg(effect_type)
    if (!identical(delimiter, "\t") && !identical(delimiter, ",")) stop("delimiter must be tab or comma")
    batch_bytes <- .mss_positive_integer(batch_bytes, "batch_bytes")
    max_record_bytes <- .mss_positive_integer(max_record_bytes, "max_record_bytes")
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
    reader <- gzfile(snapshot, open = "rb")
    on.exit(close(reader), add = TRUE)
    source_name <- if (delimiter == "\t") "source.tsv" else "source.csv"
    source_writer <- file(file.path(staging, source_name), "wb")
    on.exit(close(source_writer), add = TRUE)
    parser <- mss_parser_open_cpp(delimiter, max_record_bytes)
    pending <- list(); ended <- FALSE
    read_records <- function(n) {
        result <- list(); used <- 0
        repeat {
            if (length(pending)) {
                take <- min(n - length(result), length(pending))
                sizes <- vapply(pending[seq_len(take)], function(x) sum(nchar(x, type = "bytes")) + length(x), numeric(1))
                within <- which(cumsum(sizes) + used <= batch_bytes)
                take <- if (length(within)) max(within) else if (!length(result)) 1L else 0L
                if (!take) break
                result <- c(result, pending[seq_len(take)])
                used <- used + sum(sizes[seq_len(take)])
                pending <<- pending[-seq_len(take)]
                if (length(result) >= n || used >= batch_bytes) break
            } else if (ended) break else {
                bytes <- withCallingHandlers(readBin(reader, "raw", n = min(65536L, batch_bytes)),
                         warning = function(w) stop("Input compressed read failure: ", conditionMessage(w)))
                ended <<- !length(bytes)
                if (length(bytes)) writeBin(bytes, source_writer)
                pending <<- mss_parser_feed_cpp(parser, bytes, ended)
            }
        }
        result
    }
    header <- read_records(1L)
    if (length(header) != 1L) stop("Missing input header")
    header_names <- header[[1L]]
    if (any(!nzchar(header_names)) || anyDuplicated(header_names)) stop("Header names must be nonempty and unique")
    normalized_writer <- file(file.path(staging, "normalized.tsv"), "wt")
    on.exit(close(normalized_writer), add = TRUE)
    total <- 0
    counts <- numeric()
    projection_counts <- numeric()
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
        fields <- read_records(batch_size)
        if (!length(fields)) break
        if (any(lengths(fields) != length(header_names))) stop("Input row has an unexpected number of fields near source row ", total + 1)
        source_batch <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE, check.names = FALSE)
        names(source_batch) <- header_names
        result <- harmonise_sumstats(source_batch, reference, columns = columns, build = build,
                                     frequency_type = frequency_type, effect_type = effect_type,
                                     batch_size = batch_size, candidate_limit = candidate_limit,
                                     verify = FALSE, source_strand = source_strand, projection = projection)$normalized
        result$source_row <- result$source_row + total
        write_normalized(result)
        batch_counts <- table(result$status)
        for (status in names(batch_counts)) {
            if (!status %in% names(counts)) counts[[status]] <- 0
            counts[[status]] <- counts[[status]] + as.double(batch_counts[[status]])
        }
        if (!is.null(projection)) for (status in unique(result$projection_status)) {
            if (!status %in% names(projection_counts)) projection_counts[[status]] <- 0
            projection_counts[[status]] <- projection_counts[[status]] + sum(result$projection_status == status)
        }
        total <- total + length(fields)
    }
    if (first) {
        empty <- as.data.frame(stats::setNames(rep(list(character()), length(header_names)), header_names),
                               stringsAsFactors = FALSE, check.names = FALSE)
        write_normalized(harmonise_sumstats(empty, reference, columns = columns, build = build,
                                           frequency_type = frequency_type, effect_type = effect_type,
                                           candidate_limit = candidate_limit, verify = FALSE,
                                           source_strand = source_strand, projection = projection)$normalized)
    }
    if (!is.null(acquisition_receipt)) {
        if (!is.list(acquisition_receipt) || !length(acquisition_receipt) ||
            is.null(names(acquisition_receipt)) || anyDuplicated(names(acquisition_receipt)) ||
            any(!names(acquisition_receipt) %in% c("rows", "md5", "bytes"))) stop("Invalid acquisition receipt")
        expected <- list(rows = total, md5 = before, bytes = file.info(snapshot)$size)
        for (key in names(acquisition_receipt)) {
            if (length(acquisition_receipt[[key]]) != 1L || is.na(acquisition_receipt[[key]]) ||
                acquisition_receipt[[key]] != expected[[key]]) stop("Acquisition receipt mismatch: ", key)
        }
    }
    flush(source_writer); flush(normalized_writer)
    if (sum(counts) != total) stop("Source/native audit row count mismatch")
    audit <- data.frame(status = names(counts), count = as.double(counts))
    utils::write.table(audit, file.path(staging, "audit.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE)
    metadata <- c(schema_version = "1", core_version = reference$core_version,
                  implementation_version = "mss-preservation-4",
                  contig_mapping_implementation = if (length(reference$contig_aliases)) "mss-contigs-1" else "not_applied",
                  reference_build = build, reference_version = reference$version,
                  reference_path = reference$path, source_rows = sprintf("%.17g", total),
                  normalized_rows = sprintf("%.17g", total), input_checksum = before,
                  checksum_algorithm = "MD5", frequency_type = frequency_type,
                  effect_type = effect_type, projection = if (is.null(projection)) "not_attempted" else "attempted",
                  delimiter = if (delimiter == "\t") "tab" else "comma",
                  source_file = source_name, parser_version = "mss-delimited-1",
                  batch_bytes = as.character(batch_bytes), max_record_bytes = as.character(max_record_bytes),
                  acquisition_receipt = if (is.null(acquisition_receipt)) "not_supplied" else "verified",
                  projection_implementation = if (is.null(projection)) "not_attempted" else projection$receipt$implementation_version,
                  source_strand = if (!is.null(source_strand)) source_strand else
                      if ("source_strand" %in% names(columns)) paste0("column:", columns$source_strand) else "not_declared")
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
    if (!is.null(projection)) {
        if (sum(projection_counts) != total) stop("Projection audit row count mismatch")
        utils::write.table(data.frame(status = names(projection_counts), count = unname(projection_counts)),
                           file.path(staging, "projection_audit.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
        receipt <- unlist(projection$receipt)
        utils::write.table(data.frame(key = names(receipt), value = unname(receipt)),
                           file.path(staging, "projection.tsv"), sep = "\t", quote = TRUE,
                           row.names = FALSE, qmethod = "double")
    }
    if (length(reference$contig_aliases)) {
        if (!file.copy(file.path(reference$path, "contig_aliases.tsv"),
                       file.path(staging, "reference_contigs.tsv")))
            stop("Cannot preserve reference contig mapping")
    }
    filenames <- c(if (length(reference$contig_aliases)) "reference_contigs.tsv",
                   if (!is.null(projection)) c("projection.tsv", "projection_audit.tsv"), "source-artifact", source_name, "normalized.tsv", "audit.tsv",
                   "metadata.tsv", "columns.tsv", "reference_checksums.tsv")
    checksums <- tools::md5sum(file.path(staging, filenames))
    names(checksums) <- filenames
    utils::write.table(data.frame(filename = filenames, md5 = unname(checksums)),
                       file.path(staging, "checksums.tsv"), sep = "\t", quote = TRUE,
                       row.names = FALSE, qmethod = "double")
    manifest <- list(schema_version = 1L, core_version = reference$core_version,
                     implementation_version = "mss-preservation-4",
                     contig_mapping_implementation = if (length(reference$contig_aliases)) "mss-contigs-1" else "not_applied",
                     reference_build = build, reference_version = reference$version,
                     reference_checksums = reference$checksums, reference_path = reference$path,
                     source_rows = total, normalized_rows = total,
                     status_counts = counts, projection_status_counts = projection_counts, input_checksum = before,
                     checksums = checksums, checksum_algorithm = "MD5",
                     columns = columns, source_strand = source_strand, frequency_type = frequency_type,
                     effect_type = effect_type, projection = if (is.null(projection)) "not_attempted" else "attempted",
                     parser_version = "mss-delimited-1", batch_bytes = batch_bytes,
                     max_record_bytes = max_record_bytes, acquisition_receipt = acquisition_receipt,
                     projection_receipt = if (is.null(projection)) NULL else projection$receipt,
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
