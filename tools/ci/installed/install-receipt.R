# Current normal installation identity, not an installation substitute.
# JSON source/dependency provenance is independently checked by the receiver.
main <- function() {
    args <- commandArgs(TRUE)
    stopifnot(length(args) == 6L, args[[5L]] %in% c("linux-gold", "windows"), !file.exists(args[[6L]]))
    library <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
    dependencies <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
    .libPaths(c(library, dependencies, .Library), include.site = FALSE)
    source <- jsonlite::fromJSON(args[[3L]], simplifyVector = FALSE)
    closure <- jsonlite::fromJSON(args[[4L]], simplifyVector = FALSE)
    lane <- args[[5L]]
    stopifnot(identical(source$schema, "mss-checkout-v1"), isTRUE(source$normal_current_checkout),
              identical(as.character(getRversion()), closure$R), .Machine$sizeof.pointer == 8L,
              if (lane == "windows") identical(.Platform$OS.type, "windows") else identical(Sys.info()[["sysname"]], "Linux"))
    suppressPackageStartupMessages(library(MungeSumstats, lib.loc = library))
    actual <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "path"), winslash = "/", mustWork = TRUE)
    stopifnot(identical(actual, file.path(library, "MungeSumstats")), !pkgload::is_dev_package("MungeSumstats"))
    dll <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "DLLs")[["MungeSumstats"]][["path"]],
                         winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(dll, paste0(actual, "/libs/")))
    wanted <- c(names(Filter(function(row) identical(row$kind, "code"), closure$packages)),
                unlist(closure$lanes[[lane]]$reference_packages, use.names = FALSE))
    packages <- installed.packages(lib.loc = dependencies)
    stopifnot(setequal(rownames(packages), wanted))
    for (name in wanted) stopifnot(identical(packages[name, "Version"], closure$packages[[name]]$version),
        identical(normalizePath(find.package(name), winslash = "/", mustWork = TRUE), file.path(dependencies, name)))
    jsonlite::write_json(list(schema = "mss-current-install-v1", lane = lane,
        normal_current_checkout = TRUE, checkout_commit = source$commit, checkout_tree = source$tree,
        source_files = source$files, package_version = as.character(packageVersion("MungeSumstats")),
        library = library, package_path = actual, dll = dll, dll_md5 = unname(tools::md5sum(dll)),
        dependencies = as.list(stats::setNames(packages[wanted, "Version"], wanted)),
        R_version = as.character(getRversion()), R_home = R.home(), startup_home = path.expand("~"),
        startup_R_USER = Sys.getenv("R_USER"), startup_HOME = Sys.getenv("HOME"),
        development = FALSE, reference_payload_inventory_still_required = TRUE),
        args[[6L]], pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
}
main()
