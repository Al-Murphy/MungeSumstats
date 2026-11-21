# Run with the explicit library and runtime pins in the component proposal.
# Sourced by benchmark_flip_indices.R for shared helpers.
physical_path <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)

load_flip_runtime <- function(library_root, version, dll,
                              framework_root, framework_version, package_root) {
    library_root <- physical_path(library_root)
    framework_root <- physical_path(framework_root)
    package_root <- physical_path(package_root)
    same(physical_path(R.home()), framework_root, "approved physical R home")
    same(as.character(getRversion()), framework_version, "approved R version")
    framework_lib <- physical_path(file.path(framework_root, "library"))
    same(physical_path(.libPaths()), framework_lib, "initial R library search")
    same(list.files(library_root, all.files = TRUE, no.. = TRUE),
         "data.table", "private library package inventory")
    link <- file.path(library_root, "data.table")
    stopifnot(isTRUE(nzchar(Sys.readlink(link), keepNA = TRUE)))
    same(physical_path(link), package_root, "private library's physical data.table")
    .libPaths(c(library_root, framework_lib), include.site = FALSE)
    suppressPackageStartupMessages(
        library("data.table", lib.loc = library_root, character.only = TRUE)
    )
    # getS3method() needs tools for the primitive "[" generic.
    loadNamespace("tools", lib.loc = framework_lib)
    data.table::setDTthreads(1L)
    list(libraries = c(library_root, framework_lib), framework_version = framework_version,
         package = package_root, version = version, dll = physical_path(dll))
}

flip_runtime_receipt <- function(runtime, phase) {
    jit <- compiler::enableJIT(-1L)
    threads <- data.table::getDTthreads()
    roots <- physical_path(.libPaths())
    same(roots, runtime$libraries, "R library search roots")
    stopifnot(jit == 0L, threads == 1L)

    framework <- c("base", "compiler", "grDevices", "graphics", "methods", "stats", "tools", "utils")
    namespaces <- sort(loadedNamespaces())
    cat(sprintf("Runtime %s namespace roster: %s\n", phase, paste(namespaces, collapse = ", ")))
    same(namespaces, sort(c(framework, "data.table")), "R namespace inventory")
    expected_paths <- c(stats::setNames(file.path(roots[[2L]], framework), framework),
                        data.table = runtime$package)
    expected_versions <- c(stats::setNames(rep(runtime$framework_version, length(framework)),
                                           framework), data.table = runtime$version)
    paths <- vapply(namespaces, function(package) {
        if (identical(package, "base")) return(physical_path(file.path(.Library, "base")))
        physical_path(getNamespaceInfo(asNamespace(package), "path"))
    }, character(1L))
    versions <- vapply(namespaces, function(package) {
        as.character(getNamespaceVersion(package))
    }, character(1L))
    same(paths, expected_paths[namespaces], "R physical namespace paths")
    same(versions, expected_versions[namespaces], "R namespace versions")

    dlls <- getLoadedDLLs()
    dll_paths <- vapply(dlls, `[[`, character(1L), "path")
    framework_dlls <- c("grDevices", "graphics", "methods", "stats", "tools", "utils")
    expected_dlls <- c(stats::setNames(
        file.path(expected_paths[framework_dlls], "libs",
                  paste0(framework_dlls, .Platform$dynlib.ext)),
        framework_dlls
    ), data_table = runtime$dll)
    # The base entry describes the built-in loader, without a path to a package DLL.
    if ("base" %in% names(dll_paths)) same(dll_paths[["base"]], "base", "base DLL entry")
    physical_dlls <- dll_paths[names(dll_paths) != "base"]
    same(sort(names(physical_dlls)), sort(names(expected_dlls)), "R DLL inventory")
    physical_dlls <- vapply(physical_dlls, physical_path, character(1L))
    same(physical_dlls,
         vapply(expected_dlls[names(physical_dlls)], physical_path, character(1L)),
         "R physical DLL paths")
    dll_paths[names(physical_dlls)] <- physical_dlls

    defaults <- list(datatable.optimize = Inf, datatable.auto.index = TRUE,
                     datatable.use.index = TRUE, datatable.verbose = FALSE)
    active <- lapply(names(defaults), function(name) getOption(name, defaults[[name]]))
    same(active, unname(defaults), "data.table optimization and verbosity options")
    all_options <- options()
    dt_option_names <- sort(names(all_options)[startsWith(names(all_options), "datatable.")])
    dt_options <- modifyList(defaults, all_options[dt_option_names])
    cat(sprintf("Runtime %s: %s; R home=%s; JIT=%d; data.table threads=%d\n",
                phase, R.version.string, physical_path(R.home()), jit, threads))
    cat("Library roots:", paste(roots, collapse = " | "), "\n")
    cat(paste("Namespace", namespaces, versions, paths, sep = "\t"), sep = "\n")
    cat(paste("DLL", names(dll_paths), dll_paths, sep = "\t"), sep = "\n")
    cat("data.table options:\n")
    dput(dt_options)
    invisible(list(roots = roots, namespaces = paths, versions = versions,
                   dlls = dll_paths, jit = jit, threads = threads, options = dt_options))
}

source_flip <- function(path, helper_dir = NULL) {
    namespace <- new.env(parent = globalenv())
    if (!is.null(helper_dir)) {
        helpers <- c("messager.R", "get_unique_name_log_file.R", "write_sumstats.R")
        for (helper in helpers) {
            sys.source(file.path(helper_dir, helper), envir = namespace)
        }
    }
    sys.source(path, envir = namespace)
    namespace$check_allele_flip
}

same <- function(actual, expected, label) {
    if (!identical(actual, expected, num.eq = FALSE, single.NA = FALSE,
                   attrib.as.set = FALSE)) {
        differences <- utils::head(all.equal(actual, expected), 8L)
        cat(label, ": ", paste(differences, collapse = "\n"), "\n", sep = "", file = stderr())
        if (is.list(actual) && is.list(expected) &&
            "conditions" %in% names(actual) && "conditions" %in% names(expected)) {
            cat("Actual conditions:\n", file = stderr())
            dput(actual$conditions, file = stderr())
            cat("Expected conditions:\n", file = stderr())
            dput(expected$conditions, file = stderr())
        }
        stop(label, call. = FALSE)
    }
}

table_snapshot <- function(table) {
    if (is.null(table)) return(NULL)
    metadata <- attributes(table)
    # Self-reference pointers belong to each independently allocated data.table.
    # Keep the attribute name/order and pointer type; compare all other metadata.
    metadata[[".internal.selfref"]] <- typeof(metadata[[".internal.selfref"]])
    list(columns = lapply(table, data.table::copy), attributes = metadata)
}

invoke_flip <- function(flip, input, reference, settings = list(),
                        log_dir = NULL) {
    args <- modifyList(list(
        allele_flip_check = TRUE,
        allele_flip_drop = FALSE,
        allele_flip_z = TRUE,
        allele_flip_frq = TRUE,
        bi_allelic_filter = TRUE,
        flip_frq_as_biallelic = FALSE,
        imputation_ind = FALSE
    ), settings)
    flip(
        sumstats_dt = input, path = "fixture.tsv", ref_genome = "GRCh37",
        rsids = reference, allele_flip_check = args$allele_flip_check,
        allele_flip_drop = args$allele_flip_drop,
        allele_flip_z = args$allele_flip_z, allele_flip_frq = args$allele_flip_frq,
        bi_allelic_filter = args$bi_allelic_filter,
        flip_frq_as_biallelic = args$flip_frq_as_biallelic,
        imputation_ind = args$imputation_ind, log_folder_ind = !is.null(log_dir),
        check_save_out = list(log_folder = log_dir, extension = ".tsv", sep = "\t"),
        tabix_index = FALSE, nThread = 1L,
        log_files = list(alleles_dont_match_ref_gen = "existing-log.tsv"),
        standardise_headers = FALSE, mapping_file = NULL,
        dbSNP = NULL, dbSNP_tarball = NULL
    )
}

observe_flip <- function(flip, input, reference, settings = list(),
                         log_dir = NULL, retained_dir = NULL, expected_logs = character()) {
    conditions <- list()
    record <- function(type, condition) {
        conditions[[length(conditions) + 1L]] <<- list(
            type = type, class = class(condition),
            message = conditionMessage(condition), call = conditionCall(condition)
        )
    }
    result <- tryCatch(
        withCallingHandlers(
            invoke_flip(flip, input, reference, settings, log_dir),
            message = function(x) {
                record("message", x)
                invokeRestart("muffleMessage")
            },
            warning = function(x) {
                record("warning", x)
                invokeRestart("muffleWarning")
            }
        ),
        error = function(x) {
            record("error", x)
            NULL
        }
    )
    files <- list()
    if (!is.null(log_dir)) {
        log_names <- list.files(log_dir, all.files = TRUE, no.. = TRUE)
        same(log_names, expected_logs, "log directory: exact expected inventory")
        for (name in log_names) {
            path <- file.path(log_dir, name)
            size <- file.info(path)$size
            stopifnot(file_test("-f", path), identical(Sys.readlink(path), ""),
                      identical(physical_path(path), file.path(physical_path(log_dir), name)),
                      !is.na(size), size <= 65536)
            files[[name]] <- readBin(path, "raw", n = as.integer(size))
            stopifnot(length(files[[name]]) == size)
            stopifnot(file.rename(path, file.path(retained_dir, name)))
        }
    }
    list(
        result = if (!is.null(result)) list(
            names = names(result),
            sumstats = table_snapshot(result$sumstats_dt),
            reference = table_snapshot(result$rsids),
            log_files = result$log_files,
            aliases_input = identical(data.table::address(result$sumstats_dt),
                                      data.table::address(input)),
            aliases_reference = !is.null(reference) &&
                identical(data.table::address(result$rsids),
                          data.table::address(reference))
        ),
        caller_input = table_snapshot(input),
        caller_reference = table_snapshot(reference),
        conditions = conditions,
        log_bytes = files
    )
}

flip_input <- function(snp, a1, a2, ...) {
    data.table::data.table(SNP = snp, A1 = a1, A2 = a2,
                           ROW = seq_along(snp), ...)
}

flip_reference <- function(snp, allele = rep("A", length(snp))) {
    data.table::data.table(SNP = snp, ref_allele = allele, key = "SNP")
}

copy_case <- function(fixture) {
    list(input = data.table::copy(fixture$input),
         reference = data.table::copy(fixture$reference))
}

flip_cases <- function() {
    fixture <- function(input, reference = flip_reference(unique(input$SNP)),
                        settings = list(), log = NULL) {
        stopifnot(nrow(input) <= 64L, ncol(input) <= 64L,
                  is.null(reference) || nrow(reference) <= 64L)
        list(input = input, reference = reference, settings = settings, log = log)
    }
    ordinary <- flip_input(
        c("rs2", "rs1", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8"),
        c("G", rep("A", 7)), c("A", rep("G", 7)),
        BETA = c(-0.7, rep(0.5, 7)), OR = c(4, rep(2, 7)),
        LOG_ODDS = c(-0.4, rep(0.1, 7)), SIGNED_SUMSTAT = c(3, rep(1, 7)),
        Z = c(-2, rep(1, 7)), FRQ = c(0.2, rep(0.1, 7)),
        flipped = c(FALSE, FALSE, TRUE, rep(FALSE, 5))
    )
    partial <- flip_input(
        paste0("rs", c("06", "02", "04", "01", "03", "05", "07", "08", "09", "10")),
        c("G", "G", NA_character_, "C", NA_character_, "A", "AT", "T", "A", NA_character_),
        c(NA_character_, "A", "A", "T", "G", NA_character_, "A", "TA", "A", NA_character_),
        BETA = c(6, 2, 4, 1, 3, 5, 7, 8, 9, 10)
    )
    partial_ref <- flip_reference(
        c(paste0("rs0", 1:7), "rs09", "rs10"), c(rep("A", 6), NA_character_, "A", "A")
    )
    partial_frq <- data.table::copy(partial)
    data.table::set(partial_frq, j = "FRQ", value = rep(0.2, nrow(partial_frq)))
    new_flag <- data.table::copy(ordinary)
    data.table::set(new_flag, j = "flipped", value = NULL)
    missing_a2 <- flip_input("rs1", "G", "A", BETA = -0.7)
    data.table::set(missing_a2, j = "A2", value = NULL)
    mismatch <- flip_input(paste0("rs", 1:4), rep("C", 4), rep("T", 4), BETA = 1:4)
    dense <- flip_input(
        paste0("rs", 1:5), rep("G", 5), rep("A", 5),
        BETA = c(0, -0, 0.125, -1.5, NA_real_), OR = c(0, -0, 1, 2, NA_real_)
    )
    wide <- data.table::copy(ordinary)
    for (column in sprintf("PAYLOAD_%02d", seq_len(64L - ncol(wide)))) {
        data.table::set(wide, j = column, value = seq_len(nrow(wide)))
    }
    duplicates <- flip_input(
        c("rs2", "rs1", "rs2"), c("G", "A", "A"), c("A", "G", "A"),
        BETA = c(1, 2, 3)
    )
    data.table::setkeyv(duplicates, "SNP")
    data.table::setindexv(duplicates, "BETA")
    drop_log <- "alleles_dont_match_ref_gen_2.tsv"
    list(
        sparse = fixture(ordinary, settings = list(imputation_ind = TRUE)),
        zero_flip = fixture(flip_input(c("rs2", "rs1"), c("A", "A"), c("G", "C"),
                                       BETA = c(-0.7, 0.4))),
        dense = fixture(dense),
        drop_partial_logged = fixture(partial, partial_ref,
                                      list(allele_flip_drop = TRUE), log = drop_log),
        drop_logged_before_refusal = fixture(
            partial_frq, partial_ref,
            list(allele_flip_drop = TRUE, bi_allelic_filter = FALSE), log = drop_log
        ),
        keep_partial = fixture(partial, partial_ref),
        partial_without_mismatch = fixture(
            flip_input(c("rs1", "rs2"), c(NA_character_, "A"), c("G", NA_character_)),
            settings = list(allele_flip_drop = TRUE), log = character()
        ),
        dense_drop = fixture(mismatch, settings = list(allele_flip_drop = TRUE)),
        dense_keep = fixture(mismatch),
        duplicates = fixture(duplicates, flip_reference(c("rs1", "rs2", "rs2", "rs3"))),
        character_effect = fixture(flip_input(
            c("rs1", "rs2", "rs3"), c("G", "T", "A"), c("A", "A", "G"),
            BETA = c("-0.7", "bad", "unchanged")
        )),
        integer_effect = fixture(flip_input(
            "rs1", "G", "A", BETA = -3L, OR = 4L
        )),
        integer_tmp = fixture(flip_input(
            c("rs1", "rs2"), c("G", "A"), c("A", "G"), tmp = 1:2
        )),
        factor_flag = fixture(flip_input(
            c("rs1", "rs2"), c("G", "A"), c("A", "G"),
            BETA = c(-0.7, 0.4), flipped = factor(c("old", "old"))
        ), settings = list(imputation_ind = TRUE)),
        z_frq_disabled = fixture(ordinary, settings = list(
            allele_flip_z = FALSE, allele_flip_frq = FALSE, bi_allelic_filter = FALSE
        )),
        z_only = fixture(ordinary, settings = list(allele_flip_frq = FALSE)),
        frq_only = fixture(ordinary, settings = list(allele_flip_z = FALSE)),
        refused_frq = fixture(ordinary, settings = list(bi_allelic_filter = FALSE)),
        approximate_frq = fixture(ordinary, settings = list(
            bi_allelic_filter = FALSE, flip_frq_as_biallelic = TRUE
        )),
        absent_frq = fixture(dense, settings = list(bi_allelic_filter = FALSE)),
        new_imputation_flag = fixture(new_flag, settings = list(imputation_ind = TRUE)),
        disabled = fixture(ordinary, NULL, list(allele_flip_check = FALSE)),
        absent_a2 = fixture(missing_a2, NULL),
        empty = fixture(flip_input(character(), character(), character(), BETA = numeric())),
        wide = fixture(wide)
    )
}

literal_checks <- function(cases) {
    columns <- function(name) cases[[name]]$result$sumstats$columns
    first <- function(name, column) columns(name)[[column]][1L]
    # The keyed result has the unchanged rs1 first and the flipped rs2 second.
    sparse <- columns("sparse")
    same(sparse$SNP[1:2], c("rs1", "rs2"), "sparse: key order")
    same(sparse$A1[1:2], c("A", "A"), "sparse: A1")
    same(sparse$A2[1:2], c("G", "G"), "sparse: A2")
    same(sparse$BETA[1:2], c(0.5, 0.7), "sparse: BETA")
    same(sparse$OR[1:2], c(2, 0.25), "sparse: OR")
    same(sparse$LOG_ODDS[1:2], c(0.1, 0.4), "sparse: LOG_ODDS")
    same(sparse$SIGNED_SUMSTAT[1:2], c(1, -3), "sparse: SIGNED_SUMSTAT")
    same(sparse$Z[1:2], c(1, 2), "sparse: Z")
    same(sparse$FRQ[1:2], c(0.1, 0.8), "sparse: FRQ")
    same(sparse$flipped, c(FALSE, TRUE, TRUE, rep(FALSE, 5)), "sparse: flipped")
    same(cases$sparse$result$sumstats$attributes$sorted, "SNP", "sparse: key")
    same(cases$sparse$result$aliases_input, TRUE, "sparse: by-reference input")
    same(columns("zero_flip")$BETA, c(0.4, -0.7), "zero flip: BETA unchanged")
    same(1 / columns("dense")$BETA[1:2], c(-Inf, Inf), "dense: signed zero")
    same(columns("dense")$OR[1:2], c(Inf, -Inf), "dense: reciprocal signed zero")

    dropped <- columns("drop_partial_logged")
    same(dropped$SNP, c("rs02", "rs04", "rs05", "rs07", "rs08", "rs09"),
         "drop: exact keep set across NA, indels and overlap")
    same(dropped$BETA, c(-2, -4, 5, 7, 8, -9), "drop: recompute flip positions")
    same(dropped$A1, c("A", "A", "A", "AT", "T", "A"), "drop: A1")
    same(dropped$A2, c("G", NA_character_, NA_character_, "A", "TA", "A"),
         "drop: A2")
    same(cases$drop_partial_logged$result$aliases_input, FALSE, "drop: distinct table")
    original <- cases$drop_partial_logged$caller_input$columns
    same(original$BETA, as.double(1:10), "drop: caller kept pre-flip effects")
    same(original$match_type,
         c(NA, FALSE, NA, FALSE, TRUE, NA, TRUE, TRUE, FALSE, NA),
         "drop: caller kept temporary match classifications")
    log_name <- "alleles_dont_match_ref_gen_2.tsv"
    same(names(cases$drop_partial_logged$log_bytes), log_name, "drop: log file name")
    log_files <- cases$drop_partial_logged$result$log_files
    same(names(log_files), c("alleles_dont_match_ref_gen", "alleles_dont_match_ref_gen_2"),
         "drop: preserve existing log entry and append new name")
    same(log_files[[1L]], "existing-log.tsv", "drop: preserve existing log path")
    same(basename(log_files[[2L]]), log_name, "drop: returned new log path")
    log_text <- rawToChar(cases$drop_partial_logged$log_bytes[[log_name]])
    same(gsub("\r\n", "\n", log_text, fixed = TRUE), paste0(
        "SNP\tA1\tA2\tROW\tBETA\tref_gen_allele\tmatch_type\n",
        "rs01\tC\tT\t4\t1\tA\tNA\n"
    ), "drop: log TSV contents")
    same(columns("keep_partial")$SNP, paste0("rs", sprintf("%02d", 1:10)),
         "keep: retain mismatch and ambiguous rows")
    same(columns("keep_partial")$BETA, c(1, -2, 3, -4, 5, 6, 7, 8, -9, 10),
         "keep: overlap follows A2")
    same(columns("partial_without_mismatch")$SNP, c("rs1", "rs2"),
         "drop guard: NA is not sufficient to take drop branch")
    same(cases$partial_without_mismatch$log_bytes, list(), "drop guard: no file")
    same(columns("dense_drop")$SNP, character(), "dense mismatch: all dropped")
    same(columns("dense_keep")$BETA, 1:4, "dense mismatch: unchanged effects")

    same(columns("duplicates")$ROW, c(2L, 1L, 3L), "duplicates: key order")
    same(columns("duplicates")$BETA, c(2, -1, -3), "duplicates: update both rows")
    same(cases$duplicates$result$reference$columns$SNP, c("rs1", "rs2", "rs2"),
         "duplicates: returned reference multiplicity")
    same(columns("character_effect")$BETA, c("0.7", NA_character_, "unchanged"),
         "character effect: coercion and NA")
    stopifnot(any(vapply(cases$character_effect$conditions, function(condition) {
        condition$type == "warning" && grepl("NAs introduced", condition$message, fixed = TRUE)
    }, logical(1L))))
    same(columns("integer_effect")$BETA, 3L, "integer effect: signed BETA")
    same(columns("integer_effect")$OR, 0L, "integer effect: reciprocal assignment")

    tmp_case <- cases$integer_tmp
    stopifnot(!is.null(tmp_case$result), !"tmp" %in% names(tmp_case$caller_input$columns))
    same(tmp_case$caller_input$columns$A2, c("G", "G"), "integer tmp: caller swap")
    stopifnot(any(vapply(tmp_case$conditions, function(condition) {
        condition$type == "warning" &&
            grepl("Coercing 'character' RHS to 'integer'", condition$message, fixed = TRUE) &&
            grepl("named 'tmp'", condition$message, fixed = TRUE)
    }, logical(1L))))

    flag_case <- cases$factor_flag
    same(flag_case$result, NULL, "factor flag: no result")
    condition <- tail(flag_case$conditions, 1L)[[1L]]
    stopifnot(condition$type == "error",
              grepl("Cannot assign 'logical' to 'factor'", condition$message, fixed = TRUE))
    same(flag_case$caller_input$columns$BETA, c(0.7, 0.4), "factor flag: effects precede error")
    same(flag_case$caller_input$columns$flipped, factor(c("old", "old")),
         "factor flag: caller retains original levels and values")

    same(columns("z_frq_disabled")$Z[2L], -2, "disabled Z: unchanged")
    same(columns("z_frq_disabled")$FRQ[2L], 0.2, "disabled FRQ: unchanged")
    same(columns("z_frq_disabled")$flipped[2L], FALSE, "disabled imputation flag")
    same(columns("z_only")$Z[2L], 2, "Z only: flipped")
    same(columns("z_only")$FRQ[2L], 0.2, "Z only: FRQ unchanged")
    same(columns("frq_only")$Z[2L], -2, "FRQ only: Z unchanged")
    same(columns("frq_only")$FRQ[2L], 0.8, "FRQ only: complemented")
    same(cases$refused_frq$result, NULL, "refused FRQ: no result")
    error <- tail(cases$refused_frq$conditions, 1L)[[1L]]
    stopifnot(error$type == "error",
              grepl("only bi-allelic SNPs", error$message, fixed = TRUE))
    refused <- cases$refused_frq$caller_input$columns
    same(refused$A1[2L], "A", "refused FRQ: swap precedes the error")
    same(refused$BETA[2L], -0.7, "refused FRQ: effects follow the error")
    failed_after_log <- cases$drop_logged_before_refusal
    same(failed_after_log$result, NULL, "drop and refuse FRQ: no result")
    same(names(failed_after_log$log_bytes), log_name, "drop log precedes FRQ error")
    same(tail(failed_after_log$conditions, 1L)[[1L]]$type, "error",
         "drop and refuse FRQ: final error")
    same(failed_after_log$caller_input$columns$A1[2L], "G",
         "drop and refuse FRQ: post-drop swap leaves original alleles intact")
    same(columns("approximate_frq")$FRQ[2L], 0.8, "approximate FRQ: complemented")
    note <- tail(cases$approximate_frq$conditions, 1L)[[1L]]
    stopifnot(note$type == "message",
              grepl("non-bi-allelic SNPs", note$message, fixed = TRUE))
    same(first("absent_frq", "OR"), Inf, "absent FRQ: no refusal")
    same(columns("new_imputation_flag")$flipped, c(NA, TRUE, rep(NA, 6)),
         "new imputation flag: unflipped rows remain missing")
    same(cases$disabled$result$aliases_input, TRUE, "disabled: aliases caller")
    same(cases$disabled$conditions, list(), "disabled: silent")
    same(cases$absent_a2$result$aliases_input, TRUE, "missing A2: aliases caller")
    same(cases$absent_a2$conditions, list(), "missing A2: silent")
    same(columns("empty")$SNP, character(), "empty: typed zero rows")
    same(length(columns("wide")), 64L, "wide: all payload columns returned")
}

count_payload_subsets <- function(flip, fixture) {
    dt_namespace <- asNamespace("data.table")
    namespace_method <- function() get("[.data.table", dt_namespace, inherits = FALSE)
    original_binding <- namespace_method()
    original_dispatch <- getS3method("[", "data.table")
    stopifnot(identical(original_binding, original_dispatch),
              !inherits(original_binding, "functionWithTrace"))
    observer <- new.env(parent = emptyenv())
    observer$payload <- data.table::copy(names(fixture$input))
    observer$count <- 0L
    previous <- options(mss.flip.subsets = observer)
    on.exit({
        tryCatch({
            if (!identical(namespace_method(), original_binding)) {
                suppressMessages(untrace("[.data.table", where = dt_namespace))
            }
        }, finally = {
            tryCatch(
                registerS3method("[", "data.table", original_dispatch, envir = dt_namespace),
                finally = {
                    options(previous)
                    restored <- c(
                        binding = identical(namespace_method(), original_binding),
                        dispatch = identical(getS3method("[", "data.table"), original_dispatch),
                        option = identical(getOption("mss.flip.subsets"), previous[[1L]])
                    )
                    failed <- paste(names(restored)[!restored], collapse = ", ")
                    same(restored, c(binding = TRUE, dispatch = TRUE, option = TRUE),
                         paste("observer: incomplete restoration", failed))
                }
            )
        })
    }, add = TRUE)
    suppressMessages(trace(
        "[.data.table", where = dt_namespace, print = FALSE,
        tracer = quote({
            .mss_flip_observer <- getOption("mss.flip.subsets")
            .mss_flip_subset <- !missing(i) && missing(j) &&
                all(.mss_flip_observer$payload %in% names(x))
        }),
        exit = quote({
            if (.mss_flip_subset) {
                .mss_flip_value <- returnValue()
                if (inherits(.mss_flip_value, "data.table") &&
                    identical(names(.mss_flip_value), names(x))) {
                    .mss_flip_observer$count <- .mss_flip_observer$count + 1L
                }
            }
        })
    ))
    # The generic can dispatch through its S3 registry rather than looking up the
    # namespace binding. Register exactly the traced method for this scoped probe.
    registerS3method("[", "data.table", namespace_method(), envir = dt_namespace)

    # Calibrate actual S3 dispatch and exclude an explicit-j or narrow reference subset.
    probe <- copy_case(fixture)
    positive <- probe$input[ROW == 1L, ]
    same(nrow(positive), 1L, "observer: positive probe row count")
    same(names(positive), names(probe$input), "observer: positive probe payload")
    stopifnot(data.table::address(positive) != data.table::address(probe$input))
    positive_delta <- observer$count
    same(positive_delta, 1L, "observer: positive missing-j subset")
    invisible(probe$input[ROW == 1L, SNP])
    explicit_j_delta <- observer$count - positive_delta
    same(explicit_j_delta, 0L, "observer: explicit-j subset excluded")
    before_reference <- observer$count
    invisible(probe$reference[1L, ])
    narrow_delta <- observer$count - before_reference
    same(narrow_delta, 0L, "observer: narrow reference subset excluded")
    observer$count <- 0L
    inputs <- copy_case(fixture)
    values <- observe_flip(flip, inputs$input, inputs$reference, fixture$settings)
    list(count = observer$count, values = values,
         deltas = c(positive_delta, explicit_j_delta, narrow_delta))
}

run_flip_regression <- function(args) {
    stopifnot(length(args) == 10L)
    oracle_path <- args[[1L]]
    candidate_path <- args[[2L]]
    output <- physical_path(args[[10L]])
    stopifnot(length(list.files(output, all.files = TRUE, no.. = TRUE)) == 0L)
    runtime <- load_flip_runtime(args[[4L]], args[[5L]], args[[6L]],
                                 args[[7L]], args[[8L]], args[[9L]])
    before <- flip_runtime_receipt(runtime, "start")
    on.exit(same(flip_runtime_receipt(runtime, "end"), before,
                 "runtime closure changed during regression"), add = TRUE)
    old <- source_flip(oracle_path, args[[3L]])
    candidate <- source_flip(candidate_path, args[[3L]])
    same(formals(candidate), formals(old), "public helper formals")
    fixtures <- flip_cases()
    stopifnot(length(fixtures) <= 25L)
    outcomes <- list(control = list(), candidate = list())
    for (name in names(fixtures)) {
        fixture <- fixtures[[name]]
        log_dir <- NULL
        if (!is.null(fixture$log)) {
            log_dir <- file.path(output, name, "shared")
            stopifnot(dir.create(log_dir, recursive = TRUE, mode = "0700"))
        }
        for (role in names(outcomes)) {
            inputs <- copy_case(fixture)
            retained <- NULL
            if (!is.null(fixture$log)) {
                retained <- file.path(output, name, role)
                stopifnot(dir.create(retained, mode = "0700"))
            }
            flip <- if (role == "control") old else candidate
            outcomes[[role]][[name]] <- observe_flip(
                flip, inputs$input, inputs$reference, fixture$settings,
                log_dir, retained, expected_logs = fixture$log
            )
        }
        same(outcomes$candidate[[name]], outcomes$control[[name]], paste(name, "oracle parity"))
    }
    literal_checks(outcomes$control)
    literal_checks(outcomes$candidate)
    cat("Exact oracle and literal parity passed for", length(fixtures), "cases.\n")

    counted <- count_payload_subsets(candidate, fixtures$wide)
    same(counted$values, outcomes$candidate$wide, "observer: preserved behavior")
    cat(sprintf(paste0("Observer: missing-j=%d explicit-j=%d narrow=%d ",
                       "namespace/S3/options restored; traced parity passed.\n"),
                counted$deltas[1L], counted$deltas[2L], counted$deltas[3L]))
    same(counted$count, 0L, paste(
        "Full-width missing-j GWAS row subsets with logging and drop disabled:",
        counted$count, "(expected zero)"
    ))
    cat("Full-width missing-j GWAS row subsets with logging and drop disabled: 0\n")
}

if (sys.nframe() == 0L) {
    run_flip_regression(commandArgs(trailingOnly = TRUE))
}
