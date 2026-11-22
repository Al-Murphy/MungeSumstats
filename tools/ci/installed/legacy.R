# Run one operation per fresh, externally supervised native R process.
# Source preparation only: this file has not been parsed or run locally.
main <- function() {
    args <- commandArgs(TRUE)
    stopifnot(length(args) == 7L)
    operation <- args[[1L]]
    host <- args[[2L]]
    lib <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
    deps <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
    packet <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
    output <- args[[6L]]
    installation_path <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
    .libPaths(c(lib, deps, .Library), include.site = FALSE)
    installation <- jsonlite::fromJSON(installation_path, simplifyVector = FALSE)
    expected_dll_md5 <- installation$dll_md5
    stopifnot(operation %in% c("allele-flip", "indels", "multiple-rsids", "tabix", "headers"),
        host %in% c("windows", "linux-gold"), .Machine$sizeof.pointer == 8L,
        as.character(getRversion()) == "4.5.3",
        if (host == "windows") .Platform$OS.type == "windows" else Sys.info()[["sysname"]] == "Linux",
        grepl("^[a-f0-9]{32}$", expected_dll_md5), !file.exists(output))
    stopifnot(dir.create(output, recursive = FALSE, mode = "0700"))
    output <- normalizePath(output, winslash = "/", mustWork = TRUE)
    old <- setwd(output); on.exit(setwd(old), add = TRUE)
    .libPaths(c(lib, deps, .Library), include.site = FALSE)
    suppressPackageStartupMessages(library(MungeSumstats, lib.loc = lib))
    suppressPackageStartupMessages(library(data.table, lib.loc = deps))
    data.table::setDTthreads(1L)
    stopifnot(identical(installation$schema, "mss-current-install-v1"),
        identical(installation$lane, host), identical(installation$normal_current_checkout, TRUE),
        grepl("^[a-f0-9]{40}$", installation$checkout_commit),
        identical(as.character(packageVersion("MungeSumstats")), installation$package_version),
        normalizePath(find.package("MungeSumstats"), winslash = "/") == file.path(lib, "MungeSumstats"),
        data.table::getDTthreads() == 1L)
    dll <- normalizePath(getNamespaceInfo(getNamespace("MungeSumstats"), "DLLs")[["MungeSumstats"]][["path"]], winslash = "/", mustWork = TRUE)
    stopifnot(startsWith(dll, paste0(lib, "/MungeSumstats/libs/")),
        identical(unname(tools::md5sum(dll)), expected_dll_md5))
    input_pins <- jsonlite::fromJSON(file.path(packet, "INPUTS.json"), simplifyVector = FALSE)
    stopifnot(identical(input_pins$schema, "mss-legacy-fixtures-v1"))
    closure_path <- file.path(packet, "closure.json")
    closure <- jsonlite::fromJSON(closure_path, simplifyVector = FALSE)
    code <- names(Filter(function(item) identical(item$kind, "code"), closure$packages))
    references <- unlist(closure$lanes[[host]]$reference_packages, use.names = FALSE)
    wanted <- c(code, references)
    stopifnot(length(code) == 148L, length(wanted) == if (host == "windows") 154L else 150L,
        identical(sort(unlist(closure$legacy_used_references, use.names = FALSE)),
            sort(c("BSgenome.Hsapiens.1000genomes.hs37d5", "SNPlocs.Hsapiens.dbSNP144.GRCh37"))),
        !"seqminer" %in% wanted, setequal(names(installation$dependencies), wanted))
    installed <- installed.packages(lib.loc = deps)
    stopifnot(setequal(rownames(installed), wanted))
    for (name in wanted) stopifnot(identical(installed[name, "Version"], closure$packages[[name]]$version),
        identical(installation$dependencies[[name]], closure$packages[[name]]$version))
    # Reference source/payload hashes are independently checked by the receiver.
    # Linux's150 and Windows'154 library rosters deliberately remain distinct.
    for (name in c("BSgenome.Hsapiens.1000genomes.hs37d5", "SNPlocs.Hsapiens.dbSNP144.GRCh37",
                   "Rsamtools", "Rhtslib", "VariantAnnotation")) {
        stopifnot(normalizePath(find.package(name), winslash = "/") == file.path(deps, name))
    }
    protected <- list()
    for (name in names(input_pins$fixtures)) {
        path <- file.path(packet, name)
        pin <- input_pins$fixtures[[name]]
        stopifnot(file.info(path)$size == pin$bytes,
            identical(unname(tools::md5sum(path)), pin$md5))
        protected[[path]] <- pin$md5
    }
    for (name in c("eduAttainOkbay.txt", "ALSvcf.vcf", "ALSvcf.vcf.gz")) {
        path <- system.file("extdata", name, package = "MungeSumstats", mustWork = TRUE)
        pin <- input_pins$fixtures[[paste0("fixtures/", name)]]
        stopifnot(file.info(path)$size == pin$bytes,
            identical(unname(tools::md5sum(path)), pin$md5))
        protected[[path]] <- pin$md5
    }
    stopifnot(dir.create("inputs"), dir.create("decoded"))
    checks <- list()
    check <- function(name, value) {
        stopifnot(length(name) == 1L, !name %in% names(checks))
        checks[[name]] <<- isTRUE(value)
        if (!isTRUE(value)) stop(paste("Scientific check failed:", name))
    }
    decoded_schema <- list()
    canonical <- function(value, name, sort = "SNP") {
        value <- data.table::copy(data.table::as.data.table(value))
        if (length(sort) && all(sort %in% names(value))) data.table::setorderv(value, sort)
        data.table::fwrite(value, file.path("decoded", paste0(name, ".tsv")), sep = "\t",
            quote = TRUE, na = "NA", eol = "\n", nThread = 1L)
        decoded_schema[[name]] <<- data.frame(table = name, ordinal = seq_along(value),
            name = names(value), typeof = vapply(value, typeof, character(1L)),
            class = vapply(value, function(column) paste(class(column), collapse = ";"), character(1L)))
        value
    }
    equal <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 0))
    edu <- file.path(packet, "fixtures", "eduAttainOkbay.txt")
    raw_lines <- readLines(edu, warn = FALSE)
    check("frozen_education_fixture_93_records", length(raw_lines) == 94L)
    output_paths <- character()
    format_case <- function(input, name, ...) {
        path <- file.path(output, paste0(name, ".tsv"))
        result <- MungeSumstats::format_sumstats(input, ..., save_path = path, nThread = 1L,
            log_folder = file.path(output, paste0(name, "-logs")))
        output_paths <<- c(output_paths, path)
        result
    }
    read_table <- function(path) data.table::fread(path, nThread = 1L)
    observed <- new.env(parent = emptyenv())
    work <- function() {
        if (operation == "allele-flip") {
            changed <- raw_lines
            changed[[1L]] <- "MarkerName\tCHR\tPOS\tA2\tA1\tEAF\tBeta\tSE\tPval"
            input <- file.path("inputs", "flipped.tsv"); writeLines(changed, input, useBytes = TRUE)
            value <- read_table(input)
            value[, Beta := -Beta]; value[, EAF := 1 - EAF]
            data.table::fwrite(value, input, sep = "\t", eol = "\n", nThread = 1L)
            flipped_path <- format_case(input, "flipped", ref_genome = "GRCh37", on_ref_genome = FALSE,
                strand_ambig_filter = FALSE, bi_allelic_filter = TRUE, allele_flip_check = TRUE, dbSNP = 144)
            control_path <- format_case(edu, "control", ref_genome = "GRCh37", on_ref_genome = TRUE,
                infer_eff_direction = FALSE, strand_ambig_filter = FALSE, bi_allelic_filter = TRUE,
                allele_flip_check = TRUE, dbSNP = 144)
            control <- canonical(read_table(control_path), "control")
            flipped <- canonical(read_table(flipped_path), "flipped")
            check("historical_full_output_line_equality", setequal(readLines(flipped_path), readLines(control_path)))
            check("all_decoded_columns_and_values_equal", equal(flipped, control))
            check("nonempty_unique_variant_membership", nrow(control) > 0L && !anyDuplicated(control$SNP))
            observed$rows <- nrow(control)
            observed$expected_scope <- "Historical unmodified-input control; cross-host frozen decoded outputs required separately"
        } else if (operation == "indels") {
            source <- data.table::data.table(SNP = c("rs34589910", "rs12987662", "4:6364621"),
                CHR = c(4, 2, 4), BP = c(6364621, 100821548, 6364621), A1 = c("C", "A", "C"),
                A2 = c("CG", "C", "CG"), Uniq.a1a2 = c("4:6364621_C_CG", "aa", "4:6364621_C_CG"),
                EAF = c(.0945334, .3787, .0945334), BETA = c(-.00625732297153778, .027, -.00625732297153778),
                P = c(.4883341, 2.693e-24, .4883341))
            canonical(source, "input", sort = character())
            # Keep immutable fixture numerics before format_sumstats can mutate
            # its data.table input. This historical call has ambiguous A1/A2.
            known_input <- data.table::copy(source[SNP == "rs34589910"])
            missing_input <- data.table::copy(source[SNP == "4:6364621"])
            result <- format_case(source, "indels", ref_genome = "GRCh37", convert_small_p = TRUE,
                allele_flip_check = TRUE, snp_ids_are_rs_ids = TRUE, return_data = TRUE,
                on_ref_genome = TRUE, indels = TRUE, log_folder_ind = TRUE, dbSNP = 144)
            rows <- canonical(result$sumstats, "retained")
            check("known_indel_retained", "rs34589910" %in% rows$SNP)
            check("two_expected_rows", nrow(rows) == 2L && setequal(rows$SNP, c("rs34589910", "rs12987662")))
            indel <- rows[SNP == "rs34589910"]
            # Frozen Linux evidence confirms the documented global label
            # inference: input A1 becomes output A2, without changing BETA/EAF.
            # FRQ preserves this fixture's EAF for the inferred effect allele;
            # it is not a claim of minor-allele frequency or biological truth.
            # Indels themselves are excluded from reference-allele checks.
            check("indel_alleles_and_anchor_retained", nrow(indel) == 1L &&
                indel$CHR == 4 && indel$BP == 6364621 && indel$A1 == "CG" && indel$A2 == "C")
            check("indel_effect_and_frequency_preserved", identical(indel$BETA, known_input$BETA) &&
                identical(indel$FRQ, known_input$EAF))
            missing <- canonical(read_table(result$log_files$snp_missing_rs), "missing-rsid-log", sort = character())
            check("one_unidentified_indel_removed_and_logged", nrow(missing) == 1L &&
                missing$SNP == missing_input$SNP && missing$CHR == 4 && missing$BP == 6364621 &&
                missing$A1 == "CG" && missing$A2 == "C")
            check("logged_indel_effect_and_frequency_preserved", identical(missing$BETA, missing_input$BETA) &&
                identical(missing$FRQ, missing_input$EAF))
            observed$rows <- nrow(rows); observed$log_files <- result$log_files
        } else if (operation == "multiple-rsids") {
            changed <- raw_lines
            changed[[3L]] <- "rs9320913_rs1234_rs_45678\t6\t98584733\tA\tC\t0.5019\t0.024\t0.003\t2.457e-19"
            path <- file.path("inputs", "multiple.tsv"); writeLines(changed, path, useBytes = TRUE)
            args <- list(ref_genome = "GRCh37", on_ref_genome = FALSE, strand_ambig_filter = FALSE,
                bi_allelic_filter = FALSE, allele_flip_check = FALSE, imputation_ind = TRUE,
                remove_multi_rs_snp = FALSE, dbSNP = 144)
            altered_path <- do.call(format_case, c(list(input = path, name = "keep-first"), args))
            original_path <- do.call(format_case, c(list(input = edu, name = "single-id-control"), args))
            altered <- read_table(altered_path); original <- read_table(original_path)
            check("first_identifier_and_conversion_flag", sum(altered$convert_multi_rs_SNP %in% TRUE) == 1L &&
                altered[convert_multi_rs_SNP %in% TRUE, SNP] == "rs9320913")
            canonical(altered, "with-conversion-flag")
            altered[, convert_multi_rs_SNP := NULL]
            altered <- canonical(altered, "keep-first"); original <- canonical(original, "single-id-control")
            check("historical_full_decoded_equality", equal(altered, original))
            changed[[2L]] <- "rs12987662\t2\t100821548\tA\tC\t\t0.027\t0.003\t2.693e-24"
            changed[[4L]] <- "rs11712056\t3\t49914397\tT\tC\t0.5504\t0\t0.003\t3.304e-19"
            changed[[9L]] <- "rs2456973\t12\t56416928\tA\tC\t0.6791\t-0.02\t-0.003\t1.064e-12"
            changed[[10L]] <- "rs165633\t12\t123767929\tA\tG\t0.2257\t0.023\t0.003\t1.258e-12"
            changed[[13L]] <- "rs11191193\t4\t140764124\tA\tG\t0.6511\t0.018\t0.003\t5.444e-11"
            changed[[15L]] <- "11210860\t1\t-1\tA\tG\t0.3694\t0.017\t0.003\t2.359e-10"
            path <- file.path("inputs", "logging.tsv"); writeLines(changed, path, useBytes = TRUE)
            logged <- format_case(path, "drop-multiple", ref_genome = "GRCh37", on_ref_genome = FALSE,
                strand_ambig_filter = TRUE, bi_allelic_filter = FALSE, allele_flip_check = FALSE,
                imputation_ind = TRUE, remove_multi_rs_snp = TRUE, log_folder_ind = TRUE, dbSNP = 144)
            check("five_historical_log_categories", length(logged$log_files) == 5L)
            counts <- vapply(names(logged$log_files), function(name) {
                rows <- canonical(read_table(logged$log_files[[name]]), paste0("log-", name), sort = character())
                nrow(rows)
            }, integer(1))
            check("historical_log_row_counts", all(ifelse(grepl("snp_strand_ambiguous", names(counts)),
                counts == 8L, counts == 1L | counts == 0L)))
            observed$rows <- nrow(altered); observed$log_counts <- as.list(counts)
            canonical(read_table(logged$sumstats), "drop-multiple")
        } else if (operation == "tabix") {
            value <- MungeSumstats::formatted_example()
            path <- file.path(output, "inputs", "sorted.tsv")
            MungeSumstats::write_sumstats(value, save_path = path)
            bgz <- MungeSumstats::index_tabular(path, remove_tmp = FALSE)
            check("bgzf_and_index_exist", endsWith(bgz, ".bgz") && file.exists(bgz) && file.exists(paste0(bgz, ".tbi")))
            header <- MungeSumstats::read_header(bgz)
            check("historical_nine_column_header", ncol(header) == 9L)
            full <- MungeSumstats::read_header(bgz, n = NULL); full[, CHR := as.character(CHR)]
            check("historical_full_read_equality", equal(full, value))
            canonical(full, "full", sort = c("CHR", "BP"))
            for (method in c("GenomicRanges", "data.table")) {
                unsorted <- data.table::copy(full); data.table::setorderv(unsorted, "BETA")
                restored <- MungeSumstats:::sort_coords(data.table::copy(unsorted), sort_method = method)
                check(paste0("sort-", method), equal(restored, full) && !identical(unsorted$SNP, full$SNP))
            }
            ranges <- GenomicRanges::GRanges(seqnames = as.character(full$CHR[c(1L, nrow(full))]),
                ranges = IRanges::IRanges(start = full$BP[c(1L, nrow(full))], width = 1L))
            query <- function(ranges, label) {
                actual <- Rsamtools::scanTabix(bgz, param = ranges)
                check(paste0(label, "-range-count"), length(actual) == length(ranges))
                for (i in seq_along(actual)) {
                    expected <- full[as.character(CHR) == as.character(GenomicRanges::seqnames(ranges)[i]) &
                        BP >= IRanges::start(ranges)[i] & BP <= IRanges::end(ranges)[i]]
                    got <- data.table::fread(text = paste(actual[[i]], collapse = "\n"), sep = "\t",
                        header = FALSE, col.names = names(full), nThread = 1L)
                    got[, CHR := as.character(CHR)]
                    check(paste0(label, "-range-", i, "-decoded-equality"), equal(got, expected))
                    canonical(got, paste0(label, "-range-", i), sort = character())
                }
            }
            query(ranges, "ordered"); query(rev(ranges), "unordered")
            observed$rows <- nrow(full); observed$queried_ranges <- 4L
            observed$seqminer <- "Absent from the explicit current code/reference closure; not tested"
        } else if (operation == "headers") {
            plain_tsv <- file.path(output, "inputs", "education.tsv")
            plain_vcf <- file.path(output, "inputs", "ALS.vcf")
            stopifnot(file.copy(edu, plain_tsv), file.copy(file.path(packet, "fixtures", "ALSvcf.vcf"), plain_vcf))
            gzip_copy <- function(source, target) {
                incoming <- file(source, "rb"); outgoing <- gzfile(target, "wb")
                on.exit(close(incoming)); on.exit(close(outgoing), add = TRUE)
                repeat { raw <- readBin(incoming, "raw", 65536L); if (!length(raw)) break; writeBin(raw, outgoing) }
                target
            }
            tsv <- list(plain = plain_tsv, gzip = gzip_copy(plain_tsv, paste0(plain_tsv, ".gz")),
                        bgzf = Rsamtools::bgzip(plain_tsv, dest = paste0(plain_tsv, ".bgz"), overwrite = FALSE))
            vcf <- list(plain = plain_vcf, gzip = gzip_copy(plain_vcf, paste0(plain_vcf, ".gz")),
                        bgzf = Rsamtools::bgzip(plain_vcf, dest = paste0(plain_vcf, ".bgz"), overwrite = FALSE))
            expected_tsv <- read_table(edu)
            expected_vcf <- readLines(file.path(packet, "fixtures", "ALSvcf.vcf"), warn = FALSE)
            check("frozen_vcf_528_lines", length(expected_vcf) == 528L)
            variants <- expected_vcf[!startsWith(expected_vcf, "##")]
            expected_rows <- data.table::fread(text = paste(variants[1:3], collapse = "\n"), nThread = 1L)
            for (compression in names(tsv)) for (skip in c(FALSE, TRUE)) {
                first <- MungeSumstats::read_header(tsv[[compression]], skip_vcf_metadata = skip)
                all <- MungeSumstats::read_header(tsv[[compression]], n = NULL, skip_vcf_metadata = skip)
                label <- paste(compression, skip, sep = "-")
                # Preserve the current API's known format-dependent default:
                # plain/gzip use two input lines; BGZF uses two data records.
                wanted <- if (compression == "bgzf") 2L else 1L
                check(paste0("tsv-first-", label), equal(first, expected_tsv[seq_len(wanted)]))
                check(paste0("tsv-all-", label), equal(all, expected_tsv) && nrow(all) == 93L)
                canonical(first, paste0("tsv-first-", label), sort = character())
                canonical(all, paste0("tsv-all-", label), sort = character())
            }
            for (compression in names(vcf)) {
                metadata <- MungeSumstats::read_header(vcf[[compression]], skip_vcf_metadata = FALSE)
                rows <- MungeSumstats::read_header(vcf[[compression]], skip_vcf_metadata = TRUE)
                check(paste0("vcf-metadata-", compression), identical(metadata, expected_vcf))
                check(paste0("vcf-rows-", compression), equal(rows, expected_rows) && identical(dim(rows), c(2L, 10L)))
                writeLines(metadata, file.path("decoded", paste0("vcf-metadata-", compression, ".txt")), useBytes = TRUE)
                canonical(rows, paste0("vcf-rows-", compression), sort = character())
            }
            observed$tsv_compression_modes <- names(tsv); observed$vcf_compression_modes <- names(vcf)
            observed$current_TSV_default_data_rows <- list(plain = 1L, gzip = 1L, bgzf = 2L)
        }
    }
    started <- proc.time()
    error <- tryCatch({ work(); NULL }, error = identity)
    closing_error <- tryCatch({
        for (path in names(protected)) check(paste0("unchanged-", basename(path), "-", length(checks)),
            identical(unname(tools::md5sum(path)), protected[[path]]))
        check("loaded_dll_unchanged", identical(unname(tools::md5sum(dll)), expected_dll_md5))
        NULL
    }, error = identity)
    if (is.null(error)) error <- closing_error
    if (length(decoded_schema)) data.table::fwrite(data.table::rbindlist(decoded_schema),
        "decoded-schema.tsv", sep = "\t", quote = TRUE, na = "NA", eol = "\n", nThread = 1L)
    timing <- as.list(unclass(proc.time() - started))
    names(timing) <- names(proc.time())
    saveRDS(sessionInfo(), "session.rds")
    writeLines(capture.output(sessionInfo()), "session.txt", useBytes = TRUE)
    files <- list.files(output, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    files <- files[!file.info(files)$isdir]
    inventory <- lapply(files, function(path) list(path = substring(path, nchar(output) + 2L),
        bytes = unname(file.info(path)$size), md5 = unname(tools::md5sum(path))))
    receipt <- list(status = if (is.null(error)) "native_legacy_operation_passed" else "failed_preserved",
        operation = operation, actual_host = Sys.info()[["sysname"]], pointer_bytes = .Machine$sizeof.pointer,
        R_version = as.character(getRversion()), declared_source_commit = installation$checkout_commit,
        install_receipt_md5 = unname(tools::md5sum(installation_path)),
        dependency_code_count = length(code), installed_dependency_count = length(wanted),
        used_references = closure$legacy_used_references,
        actual_package = as.character(packageVersion("MungeSumstats")), package_path = find.package("MungeSumstats"),
        dll = dll, dll_md5 = unname(tools::md5sum(dll)), one_thread = data.table::getDTthreads() == 1L,
        dependency_closure_md5 = unname(tools::md5sum(closure_path)),
        seqminer_installed_in_qualified_library = FALSE, checks = checks, observed = as.list(observed),
        timings = timing, files = inventory,
        error = if (is.null(error)) NULL else list(class = class(error), message = conditionMessage(error)),
        closing_error = if (is.null(closing_error)) NULL else list(class = class(closing_error), message = conditionMessage(closing_error)),
        installed_source_provenance_requires_external_install_receipt = TRUE,
        cross_host_gold_comparison_still_required = TRUE)
    jsonlite::write_json(receipt, "operation.json", pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
    if (!is.null(error)) stop(error)
    cat("NATIVE_LEGACY_OPERATION_PASS", operation, "\n")
}
main()
