# Fixture setup only; all scientific expectations stay in their original cases.
mss_test_is_development <- function() {
    isNamespaceLoaded("pkgload") && isTRUE(pkgload::is_dev_package("MungeSumstats"))
}

mss_test_cleanup <- function(scratch) {
    artifacts <- Sys.getenv("MSS_TEST_ARTIFACTS", unset = "")
    if (nzchar(artifacts)) {
        artifacts <- normalizePath(artifacts, winslash = "/", mustWork = TRUE)
        target <- file.path(artifacts, basename(scratch))
        stopifnot(dir.exists(artifacts), !file.exists(target), file.copy(scratch, artifacts, recursive = TRUE))
        inventory <- function(root) {
            files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
            hashes <- tools::md5sum(files)
            names(hashes) <- substring(files, nchar(root) + 2L)
            hashes[order(names(hashes))]
        }
        stopifnot(identical(inventory(scratch), inventory(target)),
                  identical(sort(list.dirs(scratch, full.names = FALSE)), sort(list.dirs(target, full.names = FALSE))))
    }
    # Do not remove the original if evidence copying/verification failed.
    unlink(scratch, recursive = TRUE, expand = FALSE)
    stopifnot(!dir.exists(scratch))
}

mss_test_scratch <- function(pattern, .local_envir = parent.frame()) {
    scratch <- tempfile(pattern, tmpdir = tempdir())
    stopifnot(dir.create(scratch, mode = "0700"))
    withr::defer(mss_test_cleanup(scratch), envir = .local_envir)
    scratch
}

mss_test_als_fixture <- function(purpose, .local_envir = parent.frame()) {
    stopifnot(purpose %in% c("imputation", "roundtrip"))
    if (mss_test_is_development()) {
        shipped <- file.path(getNamespaceInfo(getNamespace("MungeSumstats"), "path"), "inst", "extdata", "ALSvcf.vcf")
    } else {
        shipped <- system.file("extdata", "ALSvcf.vcf", package = "MungeSumstats", mustWork = TRUE)
    }
    scratch <- mss_test_scratch(paste0("alsvcf-", purpose, "-"), .local_envir = .local_envir)
    target <- file.path(scratch, basename(shipped))
    stopifnot(file.info(shipped)$size == 47124, file.copy(shipped, target),
              identical(readBin(shipped, "raw", 47125L), readBin(target, "raw", 47125L)))
    target
}

mss_test_home_expansion <- function(.local_envir = parent.frame()) {
    scratch <- mss_test_scratch("home-expansion-", .local_envir = .local_envir)
    home <- file.path(scratch, "home"); stopifnot(dir.create(home, mode = "0700"))
    home <- normalizePath(home, winslash = "/", mustWork = TRUE)
    development <- mss_test_is_development()
    package_path <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "path"), winslash = "/", mustWork = TRUE)
    dll_path <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "DLLs")[["MungeSumstats"]][["path"]],
                              winslash = "/", mustWork = TRUE)
    dll_md5 <- unname(tools::md5sum(dll_path))
    source_dll <- if (development) file.path(package_path, "src", paste0("MungeSumstats", .Platform$dynlib.ext)) else NULL
    if (development) stopifnot(file.exists(source_dll), identical(unname(tools::md5sum(source_dll)), dll_md5))
    script <- file.path(scratch, "check-home.R"); config <- file.path(scratch, "input.rds")
    result_path <- file.path(scratch, "result.rds")
    saveRDS(list(home = home, package_path = package_path, development = development,
                 source_dll = source_dll, dll_path = dll_path, dll_md5 = dll_md5,
                 libraries = .libPaths(), result_path = result_path), config)
    writeLines(c(
        "args <- commandArgs(TRUE); stopifnot(length(args) == 1L)",
        "config <- readRDS(args[[1L]])",
        ".libPaths(config$libraries, include.site = FALSE)",
        "if (isTRUE(config$development)) {",
        "  stopifnot(requireNamespace('pkgload', quietly = TRUE), file.exists(config$source_dll), identical(unname(tools::md5sum(config$source_dll)), config$dll_md5))",
        "  pkgload::load_all(config$package_path, compile = FALSE, helpers = FALSE, attach_testthat = FALSE, export_all = FALSE, quiet = TRUE)",
        "  stopifnot(pkgload::is_dev_package('MungeSumstats'))",
        "} else {",
        "  suppressPackageStartupMessages(library(MungeSumstats, lib.loc = dirname(config$package_path)))",
        "  stopifnot(!isNamespaceLoaded('pkgload') || !pkgload::is_dev_package('MungeSumstats'))",
        "}",
        "actual <- normalizePath(getNamespaceInfo(getNamespace('MungeSumstats'), 'path'), winslash = '/', mustWork = TRUE)",
        "stopifnot(identical(actual, config$package_path))",
        "dll <- normalizePath(getNamespaceInfo(getNamespace('MungeSumstats'), 'DLLs')[['MungeSumstats']][['path']], winslash = '/', mustWork = TRUE)",
        "stopifnot(identical(unname(tools::md5sum(dll)), config$dll_md5))",
        "if (!isTRUE(config$development)) stopifnot(identical(dll, config$dll_path))",
        "stopifnot(identical(normalizePath(path.expand('~'), winslash = '/', mustWork = TRUE), config$home))",
        "out <- MungeSumstats:::check_save_path(save_path = '~/x_mss_test/a.tsv.gz', log_folder = tempdir(), log_folder_ind = FALSE, tabix_index = FALSE)",
        "saveRDS(list(home = config$home, result = out, development = config$development, actual_package_path = actual, actual_dll_path = dll, actual_dll_md5 = unname(tools::md5sum(dll)), actual_R_home = R.home()), config$result_path)"
    ), script, useBytes = TRUE)
    executable <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    candidates <- unique(c(file.path(R.home("bin"), executable), file.path(R.home(), "bin", executable)))
    candidates <- candidates[file.exists(candidates)]; stopifnot(length(candidates) > 0L)
    # Set both values before the child starts; changing HOME after startup is
    # insufficient on Windows. Restore the parent's exact set/unset state.
    old <- Sys.getenv(c("HOME", "R_USER"), unset = NA_character_, names = TRUE)
    on.exit({
        Sys.unsetenv(names(old)[is.na(old)])
        if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
    }, add = TRUE)
    Sys.setenv(HOME = home, R_USER = home)
    status <- system2(candidates[[1L]], args = c("--vanilla", shQuote(script), shQuote(config)),
                      stdout = file.path(scratch, "stdout.txt"), stderr = file.path(scratch, "stderr.txt"),
                      wait = TRUE, timeout = 60)
    if (!identical(status, 0L) || !file.exists(result_path))
        stop("Fresh R home-expansion check failed; MSS_TEST_ARTIFACTS retains evidence when configured: ", scratch)
    checked <- readRDS(result_path)
    stopifnot(identical(checked$home, home), identical(checked$actual_package_path, package_path),
              identical(checked$development, development), identical(checked$actual_dll_md5, dll_md5))
    checked
}
