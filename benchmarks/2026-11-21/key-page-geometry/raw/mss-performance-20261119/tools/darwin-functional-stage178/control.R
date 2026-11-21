args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
mode <- args[[1L]]
lib <- args[[2L]]
overlay <- args[[3L]]
.libPaths(c(lib, overlay), include.site = FALSE)
expected <- normalizePath(c(lib, overlay, file.path(R.home(), "library")), mustWork = TRUE)
stopifnot(identical(.libPaths(), expected))
options(warn = 1)

if (identical(mode, "suite")) {
    stopifnot(length(args) == 6L)
    Sys.setenv(NOT_CRAN = "true")
    source(args[[4L]], local = TRUE)
    qualify_suite(args[[5L]], args[[6L]], package = "MungeSumstats")
} else {
    stopifnot(identical(mode, "cold"), length(args) == 12L)
    stopifnot(identical(R.version.string, "R version 4.5.3 (2026-03-11)"),
              identical(R.version$platform, "aarch64-apple-darwin20"))
    first <- args[[4L]]
    out <- args[[5L]]
    expected_dll <- args[[6L]]
    references <- args[7L:12L]
    packages <- c("MungeSumstats", "testthat")
    stopifnot(first %in% packages, !any(c(packages, references) %in% loadedNamespaces()),
              !file.exists(out), dir.create(out, mode = "0700"))

    capture_state <- function(prefix) {
        namespaces <- sort(loadedNamespaces())
        metadata <- lapply(namespaces, function(name) {
            ns <- getNamespace(name)
            path <- if (identical(name, "base")) file.path(.Library, "base") else getNamespaceInfo(ns, "path")
            c(name = name, version = as.character(getNamespaceVersion(ns)),
              path = normalizePath(path, mustWork = TRUE))
        })
        libraries <- getLoadedDLLs()
        dlls <- lapply(sort(names(libraries)), function(name) {
            raw <- libraries[[name]][["path"]]
            path <- if (identical(name, "base") && identical(raw, "base")) raw else
                normalizePath(raw, mustWork = TRUE)
            c(name = name, path = path)
        })
        write_tsv(do.call(rbind, metadata), paste0(prefix, "-namespaces.tsv"))
        write_tsv(do.call(rbind, dlls), paste0(prefix, "-dlls.tsv"))
    }

    write_tsv <- function(rows, name) {
        utils::write.table(rows, file.path(out, name), sep = "\t", quote = TRUE,
                           qmethod = "double", row.names = FALSE, fileEncoding = "UTF-8")
    }

    for (name in c(first, setdiff(packages, first))) {
        message("Loading namespace ", name)
        loadNamespace(name)
    }
    stopifnot(identical(getExportedValue("data.table", "getDTthreads")(), 1L))
    dll <- getLoadedDLLs()[["MungeSumstats"]]
    stopifnot(inherits(dll, "DLLInfo"), identical(dll[["dynamicLookup"]], FALSE),
              identical(normalizePath(dll[["path"]], mustWork = TRUE), expected_dll))
    capture_state("before-references")

    resolutions <- lapply(references, function(name) {
        message("Resolving reference ", name)
        resolved <- normalizePath(find.package(name), mustWork = TRUE)
        ns <- loadNamespace(name)
        path <- normalizePath(getNamespaceInfo(ns, "path"), mustWork = TRUE)
        stopifnot(identical(resolved, path))
        version <- as.character(getNamespaceVersion(ns))
        message(name, " ", version, " ", path)
        c(name = name, version = version, path = path)
    })
    write_tsv(do.call(rbind, resolutions), "references.tsv")
    capture_state("after-references")
}
