# Run by run.py in a fresh, explicitly selected installed package.
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Expected run CONFIG or manifests OUTPUT MANIFEST...")

stable_manifest <- function(manifest) {
    if (!is.list(manifest) || !"created_utc" %in% names(manifest) ||
        length(manifest$created_utc) != 1L || is.na(manifest$created_utc))
        stop("Invalid edition manifest timestamp")
    manifest[["created_utc"]] <- NULL
    if ("path" %in% names(manifest)) manifest[["path"]] <- NULL
    manifest
}

if (identical(args[[1L]], "manifests")) {
    if (length(args) < 4L) stop("Expected output and at least two manifests")
    manifests <- lapply(args[-c(1L, 2L)], readRDS)
    if (any(vapply(manifests, function(x) "path" %in% names(x), logical(1))))
        stop("Unexpected path in on-disk edition manifest")
    manifests <- lapply(manifests, stable_manifest)
    equal <- vapply(manifests, identical, logical(1), y = manifests[[1L]])
    writeLines(if (all(equal)) "pass" else "fail", args[[2L]])
    if (!all(equal)) stop("Stable RDS manifest fields differ")
    quit(save = "no")
}
if (!identical(args[[1L]], "run") || length(args) != 2L) stop("Invalid worker command")

settings <- utils::read.delim(args[[2L]], colClasses = "character", na.strings = NULL,
                              check.names = FALSE, comment.char = "", fileEncoding = "UTF-8")
stopifnot(identical(names(settings), c("key", "value")), !anyDuplicated(settings$key))
value <- stats::setNames(as.list(settings$value), settings$key)
write_tsv <- function(x, name) {
    utils::write.table(x, file.path(value$output, name), sep = "\t", quote = TRUE,
                       qmethod = "double", row.names = FALSE, fileEncoding = "UTF-8")
}
platform_path <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
metric <- function(phase, timing) {
    stats::setNames(as.list(as.numeric(timing)), paste(phase, names(timing), sep = "."))
}

start <- proc.time()
libraries <- unlist(value[startsWith(names(value), "library.")], use.names = FALSE)
# Keep the chosen package first and use the same dependency paths in both runs.
.libPaths(libraries, include.site = FALSE)
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(MungeSumstats, lib.loc = libraries[[1L]]))
data.table::setDTthreads(1L)
package_setup <- proc.time() - start

dll <- getLoadedDLLs()[["MungeSumstats"]][["path"]]
dll_md5 <- unname(tools::md5sum(dll))
if (!identical(platform_path(dll), platform_path(value$dll)) ||
    !identical(dll_md5, value$dll_md5) ||
    !identical(platform_path(getNamespaceInfo(asNamespace("MungeSumstats"), "path")),
               platform_path(file.path(libraries[[1L]], "MungeSumstats"))))
    stop("Unexpected installed package or loaded DLL")
if (data.table::getDTthreads() != 1L) stop("Expected one data.table thread")

edition <- file.path(value$output, "edition")
if (identical(value$mode, "legacy")) {
    # Isolate the prior matched-panel adapter from the preservation API.
    # Reference table acquisition stays inside format_sumstats as in the old profile.
    load_panel <- function() {
        sites <- data.table::fread(value$legacy_panel, nThread = 1L)
        map <- Biostrings::IUPAC_CODE_MAP
        pairs <- vapply(map, function(s) paste(sort(strsplit(s, "")[[1L]]), collapse = ""),
                        character(1))
        codes <- stats::setNames(names(map), pairs)
        if (any(nchar(sites$ref) != 1L | nchar(sites$alt) != 1L | sites$ref == sites$alt))
            stop("Legacy adapter requires biallelic SNPs")
        pair <- ifelse(sites$ref < sites$alt, paste0(sites$ref, sites$alt),
                       paste0(sites$alt, sites$ref))
        ambiguity <- unname(codes[pair])
        if (anyNA(ambiguity)) stop("Legacy adapter encountered unsupported alleles")
        panel <- data.table::data.table(SNP = sites$id, seqnames = sites$chr, pos = sites$pos,
                                       alleles_as_ambig = ambiguity, ref_allele = sites$ref,
                                       alt_alleles = sites$alt)
        data.table::setkeyv(panel, "SNP")
        panel
    }
    adapter <- function(snps, ...) {
        panel <- load_panel()
        result <- panel[unique(snps), nomatch = NULL]
        data.table::setkeyv(result, "SNP")
        result
    }
    namespace <- asNamespace("MungeSumstats")
    unlockBinding("load_ref_genome_data", namespace)
    assign("load_ref_genome_data", adapter, envir = namespace)
    lockBinding("load_ref_genome_data", namespace)
    reference_open <- NULL
    process <- function() MungeSumstats::format_sumstats(
        value$input, ref_genome = "GRCh37", dbSNP = 144, on_ref_genome = FALSE,
        infer_eff_direction = FALSE, check_dups = FALSE, strand_ambig_filter = FALSE,
        frq_is_maf = FALSE, sort_coordinates = FALSE, nThread = 1,
        save_path = file.path(value$output, "legacy.tsv"), return_data = TRUE)
} else {
    start <- proc.time()
    reference <- MungeSumstats::open_reference(value$reference, fasta = value$fasta)
    reference_open <- proc.time() - start
    if (!identical(reference$core_version, value$core_version))
        stop("Unexpected reference format")
    if (!identical(platform_path(reference$fasta), platform_path(value$fasta_pin_path)))
        stop("Reference did not open the pinned FASTA")
    call <- list(input = value$input, reference = reference, path = edition,
                 frequency_type = value$frequency_type, effect_type = value$effect_type,
                 delimiter = if (identical(value$delimiter, "tab")) "\t" else ",",
                 batch_size = as.integer(value$batch_size),
                 batch_bytes = as.integer(value$batch_bytes),
                 max_record_bytes = as.integer(value$max_record_bytes),
                 candidate_limit = as.integer(value$candidate_limit),
                 build = value$build, source_strand = value$source_strand)
    columns <- value[startsWith(names(value), "column.")]
    if (length(columns)) {
        names(columns) <- substring(names(columns), 8L)
        call$columns <- columns
    }
    if (!is.null(value$acquisition_md5))
        call$acquisition_receipt <- list(rows = as.double(value$rows),
                                         md5 = value$acquisition_md5,
                                         bytes = as.double(value$acquisition_bytes))
    process <- function() do.call(MungeSumstats::harmonise_file, call)
}

# Match the prior profile's fresh-process preparation for both revisions.
heap_before <- gc(reset = TRUE)
if (identical(value$profile, "cpu"))
    utils::Rprof(file.path(value$output, "Rprof.out"), interval = 0.005,
                 memory.profiling = TRUE, gc.profiling = TRUE)
if (identical(value$profile, "alloc"))
    utils::Rprofmem(file.path(value$output, "Rprofmem.out"),
                    threshold = as.integer(value$allocation_threshold))

start <- proc.time()
result <- process()
processing <- proc.time() - start

if (identical(value$profile, "cpu")) utils::Rprof(NULL)
if (identical(value$profile, "alloc")) utils::Rprofmem(NULL)
heap_after <- gc()
if (identical(value$mode, "legacy")) {
    # Use the serialized artifact for parity and report the actual return count.
    if (!is.data.frame(result)) stop("Legacy adapter expected a returned table")
    rows <- list(legacy_returned_rows = nrow(result))
} else {
    disk <- readRDS(file.path(edition, "manifest.rds"))
    if ("path" %in% names(disk) || !identical(result$created_utc, disk$created_utc) ||
        !identical(stable_manifest(result), stable_manifest(disk)) ||
        !identical(platform_path(result$path), platform_path(edition)))
        stop("Returned manifest disagrees with published edition")
    expected <- as.double(value$rows)
    if (!identical(as.double(result$source_rows), expected) ||
        !identical(as.double(result$normalized_rows), expected))
        stop("Unexpected source or normalized row count")
    rows <- list(source_rows = result$source_rows, normalized_rows = result$normalized_rows,
                 implementation = result$implementation_version)
}

namespace_names <- sort(loadedNamespaces())
namespace_paths <- vapply(namespace_names, function(name) {
    if (identical(name, "base")) {
        platform_path(file.path(.Library, name))
    } else {
        platform_path(getNamespaceInfo(asNamespace(name), "path"))
    }
}, character(1))
write_tsv(data.frame(package = namespace_names,
                     version = vapply(namespace_names, function(name) {
                         as.character(utils::packageVersion(name, lib.loc = dirname(namespace_paths[[name]])))
                     }, character(1)),
                     path = namespace_paths), "namespaces.tsv")
runtime <- list(R_version = R.version.string, R_platform = R.version$platform,
                R_home = platform_path(R.home()), locale = Sys.getlocale(),
                BLAS = unname(extSoftVersion()["BLAS"]), LAPACK = La_version(),
                threads = data.table::getDTthreads())
if (is.null(runtime$BLAS) || is.na(runtime$BLAS)) runtime$BLAS <- "unavailable"
metadata <- c(runtime, metric("package_setup", package_setup),
              if (!is.null(reference_open)) metric("reference_open", reference_open),
              metric("processing", processing), rows,
              list(dll = platform_path(dll), dll_md5 = dll_md5))
write_tsv(data.frame(key = names(metadata), value = unlist(metadata, use.names = FALSE)), "metrics.tsv")
saveRDS(list(before = heap_before, after = heap_after), file.path(value$output, "r-heap.rds"))
if (identical(value$profile, "cpu"))
    saveRDS(summaryRprof(file.path(value$output, "Rprof.out")),
            file.path(value$output, "Rprof-summary.rds"))
