# The fresh baseline installation creates all new immutable reference editions.
args <- commandArgs(TRUE)
stopifnot(length(args) == 5L)
.libPaths(c(args[[1]], args[[2]]), include.site = FALSE)
library(MungeSumstats, lib.loc = args[[1]])
root <- args[[3]]
matched <- args[[4]]
mode <- args[[5]]
if (mode == "matched") {
    old <- open_reference(file.path(matched, "matched-index"),
                          fasta = file.path(matched, "GRCh37.ensembl104.chr22.fa"))
    stopifnot(identical(old$build, "GRCh37"), identical(old$core_version, "mss-preservation-2"))
    prepare_reference(old$fasta, file.path(old$path, "source.tsv"), file.path(root, "matched_3"),
                      build = old$build, version = old$version, contig_aliases = old$contig_aliases,
                      reference_format = "mss-preservation-3")
} else if (mode == "scientific") {
    for (format in 2:3) {
        prepare_reference(file.path(root, "scientific.fa"), file.path(root, "sites.tsv"),
                          file.path(root, paste0("scientific_", format)),
                          build = "SYNTHETIC", version = "literal-scientific-v1",
                          contig_aliases = c(chrS = "S"),
                          reference_format = paste0("mss-preservation-", format))
    }
} else stop("Expected scientific or matched")
