# Ordinary installed-package suite. Fixture/home fixes live in actual test files;
# this driver copies installed tests byte-for-byte and never rewrites a test.
main <- function() {
    args <- commandArgs(TRUE)
    stopifnot(length(args) == 5L)
    package_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
    dependencies <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
    installed_tests <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
    output <- args[[4L]]
    expected_home <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
    stopifnot(!file.exists(output), dir.create(output, mode = "0700"))
    output <- normalizePath(output, winslash = "/", mustWork = TRUE)
    .libPaths(c(package_lib, dependencies, .Library), include.site = FALSE)
    suppressPackageStartupMessages(library(MungeSumstats, lib.loc = package_lib))
    actual_package <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "path"), winslash = "/", mustWork = TRUE)
    stopifnot(identical(actual_package, file.path(package_lib, "MungeSumstats")),
              identical(installed_tests, file.path(actual_package, "tests", "testthat")),
              !pkgload::is_dev_package("MungeSumstats"), .Machine$sizeof.pointer == 8L,
              identical(normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE), expected_home),
              identical(normalizePath(Sys.getenv("R_USER"), winslash = "/", mustWork = TRUE), expected_home),
              identical(normalizePath(Sys.getenv("HOME"), winslash = "/", mustWork = TRUE), expected_home))
    dll <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "DLLs")[["MungeSumstats"]][["path"]],
                         winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(dll, paste0(actual_package, "/libs/")))
    tests <- file.path(output, "tests")
    stopifnot(dir.create(tests))
    files <- list.files(installed_tests, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    info <- file.info(file.path(installed_tests, files)); files <- files[!info$isdir]
    stopifnot(length(files) > 0L, !anyDuplicated(files))
    for (name in files) {
        target <- file.path(tests, name)
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        stopifnot(file.copy(file.path(installed_tests, name), target),
                  identical(unname(tools::md5sum(file.path(installed_tests, name))), unname(tools::md5sum(target))))
    }
    file_pins <- data.frame(path = files, md5 = unname(tools::md5sum(file.path(installed_tests, files))))
    utils::write.table(file_pins, file.path(output, "installed-test-files.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
    artifact_root <- file.path(output, "fixtures"); stopifnot(dir.create(artifact_root))
    old_artifacts <- Sys.getenv("MSS_TEST_ARTIFACTS", unset = NA_character_)
    on.exit(if (is.na(old_artifacts)) Sys.unsetenv("MSS_TEST_ARTIFACTS") else Sys.setenv(MSS_TEST_ARTIFACTS = old_artifacts), add = TRUE)
    Sys.setenv(MSS_TEST_ARTIFACTS = artifact_root)
    conditions <- list()
    record <- function(file, test, type, message) {
        conditions[[length(conditions) + 1L]] <<- data.frame(file = file, test = test, type = type,
                                                          message = message, stringsAsFactors = FALSE)
    }
    write_tsv <- function(value, name) utils::write.table(value, file.path(output, name), sep = "\t",
        quote = TRUE, qmethod = "double", row.names = FALSE, fileEncoding = "UTF-8")
    result <- tryCatch(withCallingHandlers(testthat::test_dir(tests, package = "MungeSumstats",
        load_package = "installed", reporter = "summary", stop_on_warning = FALSE, stop_on_failure = FALSE),
        warning = function(w) record("", "", "suite_warning", conditionMessage(w))),
        error = function(e) { record("", "", "suite_error", conditionMessage(e)); NULL })
    counts <- data.frame(file = character(), test = character(), nb = integer(), passed = integer(),
                         skipped = logical(), failed = integer(), error = logical(), warning = integer())
    if (length(result)) {
        frame <- as.data.frame(result); counts <- frame[names(counts)]
        for (i in seq_along(result)) for (event in result[[i]]$results) {
            if (!inherits(event, "expectation_success")) record(counts$file[i], counts$test[i],
                class(event)[1L], conditionMessage(event))
        }
    }
    write_tsv(counts, "tests.tsv")
    events <- if (length(conditions)) do.call(rbind, conditions) else
        data.frame(file = character(), test = character(), type = character(), message = character())
    write_tsv(events, "conditions.tsv")
    candidate <- nrow(counts) > 0L && sum(counts$passed) > 0L &&
        !any(counts$skipped | counts$failed > 0L | counts$error) && !any(events$type == "suite_error")
    write_tsv(data.frame(candidate_suite_pass = candidate, tests = nrow(counts), passed = sum(counts$passed),
        skipped = sum(counts$skipped), failed = sum(counts$failed), errors = sum(counts$error) + sum(events$type == "suite_error"),
        test_warnings = sum(counts$warning), suite_warnings = sum(events$type == "suite_warning"),
        source_contract_verification_required = TRUE, legacy_direct_parity_required = TRUE), "qualification.tsv")
    write_tsv(data.frame(package_path = actual_package, package_version = as.character(packageVersion("MungeSumstats")),
        dll = dll, dll_md5 = unname(tools::md5sum(dll)), R_version = as.character(getRversion()),
        R_home = R.home(), home = path.expand("~"), development = pkgload::is_dev_package("MungeSumstats")), "runtime.tsv")
    if (candidate) {
        # Query the actual byte-verified private copies retained by the ordinary
        # helper. Original scratch has already been removed by deferred cleanup;
        # reconstructed original locations and actual query paths stay distinct.
        physical <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
        tmp <- physical(file.path(dirname(expected_home), "tmp"))
        scratch <- physical(tempdir())
        stopifnot(identical(dirname(scratch), tmp),
                  all(vapply(Sys.getenv(c("TMPDIR", "TEMP", "TMP")), physical, character(1L)) == tmp))
        directories <- list.files(artifact_root, pattern = "^alsvcf-", full.names = TRUE)
        stopifnot(length(directories) == 2L)
        captured <- file.path(output, "als-fixtures"); stopifnot(dir.create(captured))
        query <- GenomicRanges::GRanges("1", IRanges::IRanges(
            start = c(1L, 1000000L), end = c(999999L, 1000000L)))
        for (purpose in c("imputation", "roundtrip")) {
            selected <- directories[startsWith(basename(directories), paste0("alsvcf-", purpose, "-"))]
            stopifnot(length(selected) == 1L)
            retained <- physical(selected)
            original_dir <- file.path(scratch, basename(selected))
            stopifnot(!dir.exists(original_dir))
            input <- file.path(retained, "ALSvcf.vcf")
            observed <- vapply(c(input = input, bgzf = paste0(input, ".bgz"),
                                tabix = paste0(input, ".bgz.tbi")), physical, character(1L))
            stopifnot(all(dirname(observed) == retained), all(file.info(observed)$size > 0L))
            tabix <- Rsamtools::TabixFile(observed[["bgzf"]], index = observed[["tabix"]])
            hits <- Rsamtools::scanTabix(tabix, param = query)
            original <- readLines(observed[["input"]], warn = FALSE)
            records <- original[!startsWith(original, "#")]
            stopifnot(length(records) == 415L, length(hits) == 2L,
                      identical(unname(hits[[1L]]), records), length(hits[[2L]]) == 0L)
            copy <- file.path(captured, purpose)
            stopifnot(dir.create(copy), all(file.copy(observed, copy)))
            writeBin(charToRaw(paste0(paste(hits[[1L]], collapse = "\n"), "\n")),
                     file.path(copy, "indexed-records.txt"))
            origin_input <- file.path(original_dir, "ALSvcf.vcf")
            origin <- c(purpose = purpose, tmp = tmp, tempdir = scratch, input = origin_input,
                        bgzf = paste0(origin_input, ".bgz"), tabix = paste0(origin_input, ".bgz.tbi"))
            utils::write.table(as.data.frame(as.list(origin)), file.path(copy, "locations.tsv"),
                sep = "\t", quote = TRUE, qmethod = "double", row.names = FALSE, fileEncoding = "UTF-8")
            actual <- c(purpose = purpose, retained_directory = retained,
                        query_input = observed[["input"]], query_bgzf = observed[["bgzf"]],
                        query_tabix = observed[["tabix"]], original_absent = "TRUE",
                        reconstructed_original_paths = "TRUE", query_records = "415", empty_range_records = "0")
            utils::write.table(as.data.frame(as.list(actual)), file.path(copy, "retained-locations.tsv"),
                sep = "\t", quote = TRUE, qmethod = "double", row.names = FALSE, fileEncoding = "UTF-8")
        }
    }
    saveRDS(sessionInfo(), file.path(output, "session.rds"))
    writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"), useBytes = TRUE)
    stopifnot(identical(unname(tools::md5sum(file.path(installed_tests, files))), file_pins$md5))
    if (!candidate) stop("Incomplete/failed installed suite; all available evidence retained")
    # The external receiver still checks exact source/case/assertion coverage,
    # controls, complete scientific gates and final package/dependency trees.
    invisible(TRUE)
}
main()
