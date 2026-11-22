# Focused installed-role control. Keep the VCFs, tabix indexes and observations.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L, .Platform$OS.type == "windows",
          .Machine$sizeof.pointer == 8L, getRversion() == "4.5.3")
physical <- function(x) normalizePath(x, winslash = "/", mustWork = TRUE)
lib <- args[[1L]]
deps <- args[[2L]]
commit <- args[[3L]]
out <- args[[4L]]
.libPaths(c(lib, deps), include.site = FALSE)
stopifnot(identical(physical(.libPaths()), physical(c(lib, deps, .Library))),
          identical(physical(find.package("MungeSumstats")), physical(file.path(lib, "MungeSumstats"))),
          grepl("^[0-9a-f]{40}$", commit), dir.create(out))

write_table <- function(x, name) {
    utils::write.table(x, file.path(out, name), sep = "\t", row.names = FALSE,
                       quote = TRUE, qmethod = "double", fileEncoding = "UTF-8")
}
sites <- data.table::data.table(
    SNP = c("rs202", "rs101"), CHR = "1", BP = c(202L, 101L),
    A1 = c("C", "A"), A2 = c("T", "G"), P = c(.001, .01),
    BETA = c(.5, .25), SE = c(.25, .125)
)
sentinel <- charToRaw("caller-owned uncompressed VCF sibling\n")
read_sentinel <- function(path) readBin(path, "raw", length(sentinel) + 1L)
stages <- function() list.files(tempdir(), "^mss-indexed-vcf-", all.files = TRUE)
stopifnot(length(stages()) == 0L)
emit_vcf <- function(mode, path, index) {
    MungeSumstats::write_sumstats(
        sumstats_dt = data.table::copy(sites), save_path = path,
        save_format = if (mode == "opengwas") "opengwas" else NULL,
        ref_genome = "GRCh37", write_vcf = TRUE, tabix_index = index, return_path = TRUE
    )
}
expected_core <- c("1:101:A:G", "1:202:C:T")
query_positions <- c(101L, 202L, 303L)
queries <- GenomicRanges::GRanges("1", IRanges::IRanges(query_positions, width = 1L))
checks <- list()
for (mode in c("normal", "opengwas")) {
    for (suffix in c("vcf", "vcf.gz", "vcf.bgz")) {
        for (indexed in c(FALSE, TRUE)) {
            label <- paste(mode, suffix, as.integer(indexed), sep = "-")
            folder <- file.path(out, label)
            stopifnot(dir.create(folder))
            sibling <- file.path(folder, "ss.vcf")
            if (indexed) writeBin(sentinel, sibling)
            requested <- file.path(folder, paste0("ss.", suffix))
            result <- emit_vcf(mode, requested, indexed)
            dest <- if (indexed) file.path(folder, "ss.vcf.bgz") else requested
            stopifnot(identical(physical(result), physical(dest)), length(stages()) == 0L)
            parsed <- VariantAnnotation::readVcf(result, genome = "GRCh37")
            gr <- SummarizedExperiment::rowRanges(parsed)
            core <- paste(as.character(GenomeInfoDb::seqnames(gr)), BiocGenerics::start(gr),
                          as.character(VariantAnnotation::ref(parsed)),
                          as.character(unlist(VariantAnnotation::alt(parsed))), sep = ":")
            stopifnot(identical(unname(core), expected_core))
            data <- VariantAnnotation::geno(parsed)
            if (mode == "opengwas") {
                stopifnot(identical(as.character(VariantAnnotation::info(parsed)$RSID), c("rs101", "rs202")),
                          identical(as.numeric(data$ES), c(.25, .5)),
                          identical(as.numeric(data$LP), c(2, 3)),
                          identical(as.numeric(data$SE), c(.125, .25)))
            } else {
                stopifnot(identical(as.character(data$SNP), c("rs101", "rs202")),
                          identical(as.numeric(data$BETA), c(.25, .5)),
                          isTRUE(all.equal(as.numeric(data$P), c(.01, .001), tolerance = 1e-6)),
                          identical(as.numeric(data$SE), c(.125, .25)))
            }
            query <- ""
            if (indexed) {
                stopifnot(identical(read_sentinel(sibling), sentinel))
                hits <- Rsamtools::scanTabix(Rsamtools::TabixFile(result), param = queries)
                stopifnot(identical(unname(lengths(hits)), c(1L, 1L, 0L)))
                tabix_core <- vapply(unname(hits), function(lines) {
                    if (!length(lines)) return("empty")
                    fields <- strsplit(lines, "\t", fixed = TRUE)[[1L]]
                    paste(fields[c(1L, 2L, 4L, 5L)], collapse = ":")
                }, character(1))
                stopifnot(identical(unname(tabix_core), c(expected_core, "empty")))
                query <- paste(tabix_core, collapse = ",")
            }
            checks[[length(checks) + 1L]] <- data.frame(
                case = label, output = basename(result), decoded = paste(core, collapse = ","),
                query = query, owned_stages = length(stages()))
        }
    }
}
write_table(do.call(rbind, checks), "checks.tsv")

# The locked bgzip must fail opening a directory as its output. Compare the
# direct native error with the package error after eliminating just input paths.
blocked <- file.path(out, "blocked")
stopifnot(dir.create(blocked))
sibling <- file.path(blocked, "ss.vcf")
dest <- paste0(sibling, ".bgz")
stopifnot(dir.create(dest))
writeBin(sentinel, sibling)
writeBin(sentinel, file.path(dest, "marker"))
catch <- function(x) tryCatch({force(x); NULL}, error = identity)
direct <- catch(Rsamtools::bgzip(sibling, dest = dest, overwrite = TRUE))
wrapped <- catch(emit_vcf("normal", sibling, TRUE))
# Registered native symbol formatting may append more than one diagnostic.
# Compare the complete messages after substituting only their input-file fields.
inputs <- function(e) {
    stopifnot(inherits(e, "error"))
    message <- conditionMessage(e)
    hits <- regmatches(message, gregexpr("(?<=\\n  file: )[^\\r\\n]+", message, perl = TRUE))[[1L]]
    stopifnot(length(hits) >= 1L, length(unique(hits)) == 1L,
              grepl("opening 'dest': ", message, fixed = TRUE))
    unique(hits)
}
a <- inputs(direct)
b <- inputs(wrapped)
raw <- chartr("\\", "/", b)
standard <- function(e, path) {
    gsub(paste0("\n  file: ", path), "\n  file: [input]", conditionMessage(e), fixed = TRUE)
}
stopifnot(identical(class(direct), class(wrapped)), identical(standard(direct, a), standard(wrapped, b)),
          identical(physical(a), physical(sibling)), identical(basename(raw), "sorted.vcf"),
          startsWith(basename(dirname(raw)), "mss-indexed-vcf-"),
          identical(physical(dirname(dirname(raw))), physical(tempdir())),
          !file.exists(dirname(raw)), length(stages()) == 0L,
          identical(read_sentinel(sibling), sentinel), identical(read_sentinel(file.path(dest, "marker")), sentinel))
write_table(data.frame(commit = commit, library = physical(lib), temp = physical(tempdir()),
                       direct = conditionMessage(direct), wrapped = conditionMessage(wrapped),
                       class = paste(class(wrapped), collapse = ","), owned_stages = length(stages())), "error.tsv")
