# Run each mode in a fresh Rscript --vanilla process with explicit libraries.
# startup checks namespace loading and semantics checks the same outputs on an
# eager baseline. Give semantics/startup a new receipt path and, optionally,
# the baseline receipt to compare every recorded field and warning.

same <- function(actual, expected, label) {
    if (!identical(actual, expected)) {
        details <- paste(all.equal(expected, actual), collapse = "; ")
        stop(label, ": ", details, call. = FALSE)
    }
}

capture_warnings <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
}

preservation <- function(directory) {
    fasta <- file.path(directory, "tiny.fa")
    sites <- file.path(directory, "sites.tsv")
    writeLines(c(">1", "ACGTACGTACGT"), fasta)
    writeLines(c(
        "chr\tpos\tref\talt\tid\taliases",
        "1\t2\tC\tT\trs2\t.",
        "1\t5\tA\tG\trs5\t."
    ), sites)
    frame <- data.frame(
        CHR = "1", BP = c(5L, 2L, 5L, 0L),
        SNP = c("rs5", "rs2", "rs5", "invalid"),
        A1 = c("A", "T", "A", "A"), A2 = c("G", "C", "G", "G"),
        BETA = c(0.4, 0.2, -0.25, 1), FRQ = c(0.2, 0.25, 0.4, 0.1),
        P = c("1e-1000", "0.2", "0.3", "0.1"), keep = letters[1:4]
    )
    expected <- data.frame(
        status = c("match", "swap", "match", "invalid_position"),
        source_row = c(1, 2, 3, 4),
        pos = c(5, 2, 5, NA_real_), effect = c(0.4, -0.2, -0.25, NA_real_),
        eaf = c(0.2, 0.75, 0.4, NA_real_)
    )
    lapply(c("mss-preservation-2", "mss-preservation-3"), function(format) {
        ref <- MungeSumstats::prepare_reference(
            fasta, sites, file.path(directory, format), build = "GRCh38",
            version = "deferral-test", reference_format = format
        )
        handle <- MungeSumstats::open_reference(ref)
        results <- lapply(list(frame, data.table::as.data.table(frame)), function(input) {
            batch <- function(size) {
                out <- MungeSumstats::harmonise_sumstats(
                    input, handle, frequency_type = "EAF", batch_size = size
                )
                same(out$source, input, "source rows, classes and tokens")
                same(out$normalized[names(expected)], expected, "resident row order and effects")
                list(
                    source_class = class(out$source), normalized = out$normalized,
                    provenance = out$provenance[c("core_version", "implementation_version",
                                                  "reference_version", "build", "frequency_type",
                                                  "effect_type", "columns", "reference_verified")]
                )
            }
            split <- batch(1L)
            same(split, batch(10L), "resident one-row batches versus a single batch")
            split
        })
        same(results[[2L]]$normalized, results[[1L]]$normalized, "data.table versus data.frame")
        results
    })
}

mapping <- function() {
    lookup <- data.table::data.table(
        uncorrected = c("pos0", "chr0", "maf"), corrected = c("BP", "CHR", "FRQ")
    )
    original <- data.table::copy(lookup)
    input <- data.frame(
        pos0 = c(7L, 2L), chr0 = c("2", "1"), maf = c(0.1, 0.2),
        A2Freq = c(0.4, 0.5), unchangedCase = c("First", "Second")
    )
    expected <- data.frame(
        BP = c(7L, 2L), CHR = c("2", "1"), maf = c(0.1, 0.2),
        FRQ = c(0.4, 0.5), unchangedCase = c("First", "Second")
    )
    out <- MungeSumstats::standardise_header(
        data.table::as.data.table(input), mapping_file = lookup,
        uppercase_unmapped = FALSE, return_list = FALSE
    )
    same(class(out), c("data.table", "data.frame"), "mapped input class")
    same(as.data.frame(out), expected, "mapping precedence, column names, classes and order")
    same(lookup, original, "caller-owned mapping")
    as.data.frame(out)
}

# The packaged ALS fixture covers the standard legacy path. These two rows
# supply ANN and positions in separate tiles for the split/merge checks.
local_vcf <- function(directory) {
    path <- file.path(directory, "annotations.vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=100,assembly=GRCh37>",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Consequence | Score | Impact\'">',
        '##FORMAT=<ID=LP,Number=1,Type=Float,Description="Negative log10 p-value">',
        paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "study"),
              collapse = "\t"),
        "1\t2\trs2\tC\tT\t12\tPASS\tANN=missense|11|moderate\tLP\t2",
        "1\t90\trs90\tA\tG\t14\tPASS\tANN=synonymous|22|low\tLP\t3"
    ), path)
    path
}

vcf_rows <- function(vcf) {
    same(methods::is(vcf, "CollapsedVCF"), TRUE, "VCF S4 class")
    ranges <- MatrixGenerics::rowRanges(vcf)
    ids <- c("rs2", "rs90")
    same(names(ranges), ids, "VCF variant IDs")
    same(GenomicRanges::start(ranges), c(2L, 90L), "VCF coordinate order and type")
    same(VariantAnnotation::samples(VariantAnnotation::header(vcf)), "study", "VCF sample header")
    same(as.character(VariantAnnotation::ref(vcf)), c("C", "A"), "VCF REF")
    if (!"ANN" %in% names(VariantAnnotation::info(vcf))) {
        stop("VCF INFO metadata lost ANN", call. = FALSE)
    }
    converted <- MungeSumstats::vcf2df(vcf, add_sample_names = FALSE, verbose = FALSE)
    same(class(converted), c("data.table", "data.frame"), "VCF conversion class")
    same(converted$ID, ids, "converted VCF IDs")
    same(converted$ANN_Consequence, c("missense", "synonymous"), "ANN consequence")
    same(converted$ANN_Score, c(11L, 22L), "ANN numeric annotation")
    as.data.frame(converted)
}

legacy_vcf <- function(directory, path) {
    # No VariantAnnotation call precedes this MungeSumstats entry point.
    vcf <- MungeSumstats::read_vcf(
        path, as_datatable = FALSE, download = FALSE, use_params = FALSE, verbose = FALSE
    )
    stopifnot("VariantAnnotation" %in% loadedNamespaces())
    stopifnot(!"package:VariantAnnotation" %in% search())
    converted <- vcf_rows(vcf)
    halves <- list(vcf[1L, ], vcf[2L, ])
    same(vcf_rows(do.call(VariantAnnotation::rbind, halves)), converted, "S4 split and merge")

    indexed <- MungeSumstats:::index_vcf(path, verbose = FALSE)
    stopifnot(file.exists(indexed), file.exists(paste0(indexed, ".tbi")))
    tabular <- MungeSumstats::read_vcf(
        indexed, download = FALSE, use_params = FALSE, verbose = FALSE
    )
    same(tabular$SNP, c("rs2", "rs90"), "legacy tabular variant IDs")
    same(tabular$BP, c(2L, 90L), "legacy tabular coordinates")
    same(tabular$P, c(0.01, 0.001), "legacy p-value conversion")

    sumstats <- data.table::data.table(
        CHR = "1", BP = c(90L, 2L), SNP = c("rs90", "rs2"),
        A1 = c("A", "C"), A2 = c("G", "T"), P = c(0.001, 0.01),
        BETA = c(-0.25, 0.5), SE = c(0.125, 0.25)
    )
    written <- MungeSumstats::write_sumstats(
        sumstats, file.path(directory, "written.vcf"), write_vcf = TRUE,
        ref_genome = "GRCh37", save_format = "opengwas", tabix_index = TRUE, return_path = TRUE
    )
    stopifnot(file.exists(written), file.exists(paste0(written, ".tbi")))
    reread <- MungeSumstats::read_vcf(
        written, download = FALSE, use_params = FALSE, verbose = FALSE
    )
    same(reread$BP, c(2L, 90L), "writer coordinate order")
    same(reread$ES, c(0.5, -0.25), "writer effect metadata")
    same(reread$RSID, c("rs2", "rs90"), "writer variant metadata")
    list(annotated = converted, indexed_read = as.data.frame(tabular), written_read = as.data.frame(reread))
}

bp_ranges <- function() {
    lapply(c("GRCh37", "GRCh38"), function(build) {
        input <- data.table::data.table(
            CHR = c("chr1", "1", "23", "chrX", "unknown", "1"),
            BP = c("5", "1000000000", "7", "0", "1000000000", "non-numeric"),
            SNP = paste0("row", 1:6)
        )
        out <- MungeSumstats:::check_bp_range(
            input, path = "input.tsv", ref_genome = build, log_folder_ind = FALSE,
            imputation_ind = FALSE, check_save_out = NULL, tabix_index = FALSE,
            nThread = 1L, log_files = list()
        )
        expected <- data.frame(
            CHR = c("chr1", "23", "unknown"), BP = c(5, 7, 1000000000),
            SNP = c("row1", "row3", "row5")
        )
        same(class(out$sumstats_dt), c("data.table", "data.frame"), "BP result class")
        same(as.data.frame(out$sumstats_dt), expected, "BP validation, classes and row order")
        same(out$log_files, list(), "BP log files")
        expected
    })
}

main <- function(args) {
    if (length(args) < 3L || length(args) > 5L ||
        !args[[1L]] %in% c("startup", "semantics", "parallel-vcf", "bp")) {
        stop("Usage: Rscript --vanilla script.R startup|semantics|parallel-vcf|bp MSS_LIB DEP_LIB [RECEIPT [BASELINE_RECEIPT]]")
    }
    mode <- args[[1L]]
    physical <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
    roots <- physical(c(args[2:3], .Library))
    .libPaths(roots, include.site = FALSE)
    same(physical(.libPaths()), roots, "installed library paths")
    stopifnot(!"VariantAnnotation" %in% loadedNamespaces())
    stopifnot(!"MungeSumstats" %in% loadedNamespaces())
    stopifnot(requireNamespace("MungeSumstats", lib.loc = roots[[1L]], quietly = TRUE))
    same(
        physical(getNamespaceInfo(asNamespace("MungeSumstats"), "path")),
        physical(file.path(roots[[1L]], "MungeSumstats")), "MungeSumstats installation"
    )
    if (mode == "startup" && "VariantAnnotation" %in% loadedNamespaces()) {
        stop("VariantAnnotation loaded while loading MungeSumstats", call. = FALSE)
    }
    data.table::setDTthreads(1L)
    directory <- tempfile("mss-namespace-")
    stopifnot(dir.create(directory))
    on.exit(unlink(directory, recursive = TRUE))
    if (mode == "bp") {
        receipt <- capture_warnings(bp_ranges())
    } else if (mode == "parallel-vcf") {
        if (!requireNamespace("GenomicFiles", quietly = TRUE)) {
            stop("The production parallel VCF check requires GenomicFiles", call. = FALSE)
        }
        receipt <- capture_warnings({
            path <- MungeSumstats:::index_vcf(local_vcf(directory), verbose = FALSE)
            vcf_rows(MungeSumstats:::read_vcf_parallel(
                path, as_datatable = FALSE, use_params = FALSE, nThread = 2L, ntile = 2L,
                verbose = FALSE
            ))
        })
    } else {
        path <- local_vcf(directory)
        before <- capture_warnings(list(preservation = preservation(directory), mapping = mapping()))
        if (mode == "startup" && "VariantAnnotation" %in% loadedNamespaces()) {
            stop("VariantAnnotation loaded during tabular processing", call. = FALSE)
        }
        vcf <- capture_warnings(legacy_vcf(directory, path))
        after <- capture_warnings(list(preservation = preservation_again(directory), mapping = mapping()))
        same(after, before, "tabular outputs and warnings after first VCF use")
        again <- capture_warnings(vcf_rows(MungeSumstats::read_vcf(
            path, as_datatable = FALSE, download = FALSE, use_params = FALSE, verbose = FALSE
        )))
        same(again$value, vcf$value$annotated, "VCF conversion after preservation")
        receipt <- list(tabular = before, vcf = vcf, repeated_vcf = again)
    }

    if (length(args) >= 5L) {
        same(receipt, readRDS(args[[5L]]), "baseline outputs and warnings")
    }
    if (length(args) >= 4L) {
        stopifnot(!file.exists(args[[4L]]))
        saveRDS(receipt, args[[4L]])
    }
    cat(mode, "passed\n")
}

preservation_again <- function(directory) {
    next_directory <- file.path(directory, "after-vcf")
    stopifnot(dir.create(next_directory))
    preservation(next_directory)
}

main(commandArgs(trailingOnly = TRUE))
