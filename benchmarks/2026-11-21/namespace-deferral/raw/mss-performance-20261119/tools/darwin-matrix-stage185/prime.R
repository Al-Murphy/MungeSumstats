# Inserted into the pinned preservation worker immediately before heap reset.
if (!identical(value$mode, "matrix") || !identical(value$profile, "none") ||
    !inherits(reference, "mss_reference_handle")) stop("Expected an unprofiled preservation handle")
prime_dir <- file.path(value$output, "prime")
if (!dir.create(prime_dir, mode = "0700")) stop("Need a fresh prime directory")
prime_path <- file.path(prime_dir, "edition")
prime_call <- call
prime_call$path <- prime_path
prime_start <- proc.time()
primed <- do.call(MungeSumstats::harmonise_file, prime_call)
prime_elapsed <- proc.time() - prime_start
prime_disk <- readRDS(file.path(prime_path, "manifest.rds"))
if ("path" %in% names(prime_disk) || !identical(primed$created_utc, prime_disk$created_utc) ||
    !identical(stable_manifest(primed), stable_manifest(prime_disk)) ||
    !identical(platform_path(primed$path), platform_path(prime_path)) ||
    !identical(primed[["reference_fasta_binding"]], reference[["fasta_binding"]]) ||
    as.double(primed$source_rows) != as.double(value$rows) ||
    as.double(primed$normalized_rows) != as.double(value$rows)) stop("Incomplete full-input prime")
write_metrics(c(metric("prime", prime_elapsed), list(source_rows = sprintf("%.0f", primed$source_rows),
                normalized_rows = sprintf("%.0f", primed$normalized_rows),
                scope = "complete identical input on the same verified handle")),
              file.path(prime_dir, "metrics.tsv"))
rm(primed, prime_disk, prime_call)
# The frozen code now resets GC accounting and times the separate full call.
