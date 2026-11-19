# Execute the worker's real open block with the public R verifier and a tiny
# synthetic reference. Stand in only for the native stream constructor.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Expected WORKER.R PRESERVATION.R NEW_DIRECTORY")
worker <- parse(args[[1L]])
scratch <- args[[3L]]
if (!dir.create(scratch)) stop("Expected a fresh synthetic fixture directory")
refdir <- file.path(scratch, "reference")
dir.create(refdir)
original <- file.path(scratch, "recorded.fa")
local_fasta <- file.path(scratch, "runtime.fa")
writeLines(c(">1", "AC"), local_fasta, useBytes = TRUE)
writeLines("1\t2\t3\t2\t3", file.path(refdir, "fasta.fai"), useBytes = TRUE)
writeLines("mss-preservation-2", file.path(refdir, "core-version"), useBytes = TRUE)
files <- c("fasta.fai", "core-version")
checks <- c(fasta = local_fasta, stats::setNames(file.path(refdir, files), files))
checksums <- stats::setNames(unname(tools::md5sum(checks)), names(checks))
sizes <- stats::setNames(as.numeric(file.info(checks)$size), names(checks))
manifest <- list(schema_version = 1L, core_version = "mss-preservation-2",
                 fasta = original, files = files, sizes = sizes, checksums = checksums)
saveRDS(manifest, file.path(refdir, "manifest.rds"))

reference_env <- new.env()
sys.source(args[[2L]], envir = reference_env)
native <- function(path, fasta) {
    stopifnot(file.exists(fasta), identical(path, normalizePath(refdir)))
    "synthetic native handle"
}
reference_env$mss_open_reference_cpp <- native

is_assignment <- function(expr, name) {
    is.call(expr) && identical(expr[[1L]], quote(`<-`)) && identical(expr[[2L]], as.name(name))
}
platform_path <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
open_branch <- Filter(function(expr) {
    is.call(expr) && length(expr) == 4L && identical(expr[[1L]], quote(`if`)) &&
        identical(expr[[2L]], quote(identical(value$mode, "legacy"))) &&
        any(vapply(as.list(expr[[4L]])[-1L], is_assignment, logical(1), name = "call"))
}, as.list(worker))
stopifnot(length(open_branch) == 1L, length(open_branch[[1L]]) == 4L)
statements <- as.list(open_branch[[1L]][[4L]])[-1L]
first_file_call <- which(vapply(statements, is_assignment, logical(1), name = "call"))
stopifnot(length(first_file_call) == 1L, first_file_call > 1L)
statements <- statements[seq_len(first_file_call - 1L)]

rewrite_namespace <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1L]], quote(MungeSumstats::open_reference)))
        expr[[1L]] <- as.name("open_reference")
    for (i in seq_along(expr)) {
        if (!identical(expr[[i]], quote(expr = )))
            expr[i] <- list(rewrite_namespace(expr[[i]]))
    }
    expr
}
run_open <- function(..., pinned = local_fasta) {
    env <- new.env(parent = reference_env)
    env$value <- list(reference = refdir, core_version = "mss-preservation-2",
                      fasta_pin_path = pinned, ...)
    env$platform_path <- platform_path
    for (statement in statements) eval(rewrite_namespace(statement), env)
    env$reference
}

explicit <- run_open(fasta = local_fasta)
stopifnot(
    identical(explicit$fasta_binding$implementation_version, "mss-fasta-binding-1"),
    identical(explicit$fasta_binding$recorded_path, original),
    identical(platform_path(explicit$fasta_binding$runtime_path), platform_path(local_fasta))
)
missing_default <- tryCatch(run_open(), error = identity)
if (!file.copy(local_fasta, original)) stop("Could not restore synthetic recorded FASTA")
default <- run_open(pinned = original)
if (!inherits(missing_default, "error"))
    stop("A missing recorded FASTA was silently replaced in default mode")
stopifnot(grepl("Reference checksum/size mismatch", conditionMessage(missing_default), fixed = TRUE))
if (!is.null(default[["fasta_binding"]]))
    stop("Default-open emitted a FASTA relocation receipt")
cat("Worker honors default and explicit FASTA opening and receipts\n")

# The revised worker exposes this base-R probe. Exercise the actual observer
# inside the same direct R process as the reference regression.
probe_definition <- Filter(function(expr) is_assignment(expr, "runtime_identity"), as.list(worker))
if (length(probe_definition)) {
    stopifnot(length(probe_definition) == 1L)
    probe_env <- new.env()
    probe_env$platform_path <- platform_path
    eval(probe_definition[[1L]], probe_env)
    observed <- probe_env$runtime_identity()
    stopifnot(identical(observed$R_home, platform_path(R.home())),
              identical(observed$R_bin, platform_path(R.home("bin"))),
              file.exists(observed$R_engine), startsWith(observed$R_engine, paste0(observed$R_bin, "/")))
    for (name in c("R_shared", "BLAS_library", "LAPACK_library"))
        stopifnot(identical(observed[[name]], "unavailable") || file.exists(observed[[name]]))
    cat("Worker observes the running R engine and available library paths\n")
}
