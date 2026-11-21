# Reconstruct the installed converter inputs from the exact phase 1 VCF.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
physical <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
roots <- physical(c(args[1:2], .Library))
.libPaths(roots, include.site = FALSE)
stopifnot(identical(physical(.libPaths()), roots))
stopifnot(requireNamespace("MungeSumstats", lib.loc = roots[[1L]], quietly = TRUE))

test_definitions <- as.list(parse(file = args[[3L]]))
fixture <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], quote(`<-`)) &&
        identical(expr[[2L]], quote(local_vcf))
}, test_definitions)
stopifnot(length(fixture) == 1L)
eval(fixture[[1L]])
path <- local_vcf(tempdir())

vcf <- MungeSumstats::read_vcf(
    path, as_datatable = FALSE, download = FALSE, use_params = FALSE, verbose = FALSE
)
expanded <- VariantAnnotation::expand(vcf)
info <- VariantAnnotation::info(expanded)
header <- VariantAnnotation::info(VariantAnnotation::header(expanded))["ANN", ]$Description
cat("ANN header and imported INFO representation\n")
dput(header)
str(info[["ANN"]])

ranges <- MatrixGenerics::rowRanges(expanded)
df <- data.table::data.table(
    ID = names(ranges), MungeSumstats:::granges_to_dt(ranges), MungeSumstats:::DF_to_dt(info)
)
cat("Converted ANN column\n")
str(df$ANN)
cat("Data table dimensions and columns before ANN removal\n")
print(dim(df))
print(names(df))

anncols <- strsplit(sub("Functional annotations: '", "", header), " \\| ")[[1L]]
empty <- data.frame(matrix(vector(), 0L, length(anncols), dimnames = list(NULL, anncols)),
                    stringsAsFactors = FALSE)
bind_args <- c(empty, lapply(lapply(df$ANN, `[`, 1), function(x) strsplit(x, "\\|")[[1L]]))
cat("Exact .anncols rbind arguments\n")
str(bind_args)

describe <- function(label, value) {
    cat(label, "\n")
    tryCatch({
        print(class(value))
        print(dim(value))
        if (is.atomic(value) && is.null(dim(value))) str(value)
        invisible(value)
    }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        invisible(NULL)
    })
}

current_rbind <- get("rbind", envir = environment(MungeSumstats::vcf2df), inherits = TRUE)
cat("MSS resolves VariantAnnotation rbind:", identical(current_rbind, VariantAnnotation::rbind), "\n")
current <- describe("Current rbind result", suppressWarnings(do.call(current_rbind, bind_args)))
base <- describe("Base rbind result", suppressWarnings(do.call(base::rbind, bind_args)))
cat("Identical current and base bind results:", identical(current, base), "\n")

original <- describe("Installed converter's ANN removal", df[, colnames(df) != "ANN"])
explicit <- describe("Explicit data.table column selection", df[, colnames(df) != "ANN", with = FALSE])
if (!is.null(current)) {
    annotations <- as.data.frame(current)
    describe("Cbind after installed converter ANN removal", cbind(original, annotations))
    describe("Cbind after explicit selection", cbind(explicit, annotations))
}
