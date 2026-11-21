# Run as a separate process for each admitted role, repetition and measurement.
make_timing_fixture <- function(width, density) {
    n <- 100000L
    row <- seq_len(n)
    flip_period <- if (density == "sparse") 100L else 4L
    flip <- if (density == "sparse") row %% flip_period == 1L else row %% flip_period != 0L
    mismatch <- !flip & row %% 100L == 0L
    key <- ifelse(flip, "rs2", ifelse(mismatch, "rs3", "rs1"))
    input <- flip_input(
        key, ifelse(flip, "G", ifelse(mismatch, "C", "A")),
        ifelse(flip, "A", ifelse(mismatch, "T", "G")),
        BETA = rep(-0.7, n), OR = rep(4, n), LOG_ODDS = rep(-0.4, n),
        SIGNED_SUMSTAT = rep(3, n), Z = rep(-2, n), FRQ = rep(0.2, n),
        flipped = rep(FALSE, n)
    )
    if (width == "wide") {
        for (column in sprintf("PAYLOAD_%02d", seq_len(16L))) {
            data.table::set(input, j = column, value = as.double(row))
        }
    }
    data.table::setkeyv(input, "SNP")
    list(input = input, reference = flip_reference(c("rs1", "rs2", "rs3", "rs4")),
         settings = list(imputation_ind = TRUE))
}

check_timing_result <- function(result, input, fixture, width, density) {
    expected_flips <- if (density == "sparse") 1000L else 75000L
    n <- nrow(input)
    result_dt <- result$sumstats_dt
    same(names(result), c("sumstats_dt", "rsids", "log_files"), "timing result fields")
    same(nrow(result_dt), 100000L, "timing row count")
    same(ncol(result_dt), if (width == "wide") 27L else 11L, "timing column count")
    same(data.table::key(result_dt), "SNP", "timing key")
    same(data.table::address(result_dt), data.table::address(input), "timing caller alias")
    same(result_dt$ROW, fixture$input$ROW, "timing row order")
    same(result_dt$SNP, fixture$input$SNP, "timing SNPs")
    same(result_dt$flipped, result_dt$SNP == "rs2", "timing flipped rows")
    same(sum(result_dt$flipped), expected_flips, "timing flip count")
    same(sum(result_dt$SNP == "rs3"), 1000L, "timing retained mismatch count")
    same(result_dt$A1, ifelse(result_dt$SNP == "rs3", "C", "A"), "timing A1")
    same(result_dt$A2, ifelse(result_dt$SNP == "rs3", "T", "G"), "timing A2")
    same(result_dt$BETA, ifelse(result_dt$flipped, 0.7, -0.7), "timing BETA")
    same(result_dt$OR, ifelse(result_dt$flipped, 0.25, 4), "timing OR")
    same(result_dt$LOG_ODDS, ifelse(result_dt$flipped, 0.4, -0.4), "timing LOG_ODDS")
    same(result_dt$SIGNED_SUMSTAT, ifelse(result_dt$flipped, -3, 3), "timing SIGNED_SUMSTAT")
    same(result_dt$Z, ifelse(result_dt$flipped, 2, -2), "timing Z")
    same(result_dt$FRQ, ifelse(result_dt$flipped, 0.8, 0.2), "timing FRQ")
    same(result$rsids$SNP, c("rs1", "rs2", "rs3"), "timing filtered reference")
    same(result$log_files, list(alleles_dont_match_ref_gen = "existing-log.tsv"),
         "timing no output log")
    if (width == "wide") {
        for (column in sprintf("PAYLOAD_%02d", seq_len(16L))) {
            same(result_dt[[column]], fixture$input[[column]], "timing payload")
        }
    }
    list(rows = n, flipped = expected_flips)
}

verify_flip_cell <- function(flip, other, fixture, width, density) {
    inputs <- copy_case(fixture)
    old <- observe_flip(flip, inputs$input, inputs$reference, fixture$settings)
    inputs <- copy_case(fixture)
    updated <- observe_flip(other, inputs$input, inputs$reference, fixture$settings)
    same(updated, old, paste(width, density, "timing oracle parity"))
    stopifnot(!is.null(updated$result))
    same(sum(updated$result$sumstats$columns$flipped),
         if (density == "sparse") 1000L else 75000L, "timing verify flips")
    same(vapply(updated$conditions, `[[`, character(1L), "type"),
         c("message", "message"), "timing verify conditions")
    cat(sprintf("%s/%s exact oracle parity\n", width, density))
}

measure_flip <- function(mode, flip, fixture, path, width, density) {
    inputs <- copy_case(fixture)
    run <- function() suppressMessages(
        invoke_flip(flip, inputs$input, inputs$reference, fixture$settings)
    )
    invisible(gc())
    if (mode == "time") {
        start <- proc.time()
        result <- run()
        elapsed <- (proc.time() - start)[c("user.self", "sys.self", "elapsed")]
        metrics <- paste(sprintf("%s=%.6f", names(elapsed), elapsed), collapse = " ")
    } else {
        Rprofmem(path, threshold = 0)
        on.exit(Rprofmem(NULL), add = TRUE)
        result <- run()
        Rprofmem(NULL)
        stopifnot(file.info(path)$size < 2 * 1024^2)
        lines <- readLines(path, warn = FALSE)
        allocated <- grep("^[0-9]+ ", lines, value = TRUE)
        bytes <- as.double(sub(" .*", "", allocated))
        stopifnot(length(bytes) > 0L, all(is.finite(bytes)))
        metrics <- sprintf("allocation_events=%d allocation_bytes=%.0f",
                           length(bytes), sum(bytes))
    }
    receipt <- check_timing_result(result, inputs$input, fixture, width, density)
    cat(sprintf("%s/%s rows=%d flipped=%d %s\n", width, density,
                receipt$rows, receipt$flipped, metrics))
}

run_flip_benchmark <- function(args) {
    stopifnot(length(args) == 11L)
    mode <- match.arg(args[[1L]], c("verify", "time", "alloc"))
    runtime <- load_flip_runtime(args[[6L]], args[[7L]], args[[8L]],
                                 args[[9L]], args[[10L]], args[[11L]])
    before <- flip_runtime_receipt(runtime, "start")
    on.exit(same(flip_runtime_receipt(runtime, "end"), before,
                 "runtime closure changed during benchmark"), add = TRUE)
    if (mode == "alloc") stopifnot(capabilities("profmem"))
    output <- physical_path(args[[5L]])
    stopifnot(length(list.files(output, all.files = TRUE, no.. = TRUE)) == 0L)
    flip <- source_flip(args[[3L]])
    other <- if (mode == "verify") source_flip(args[[4L]]) else NULL
    matrix <- list(c("narrow", "sparse"), c("wide", "dense"),
                   c("narrow", "dense"), c("wide", "sparse"))
    for (cell in matrix) {
        width <- cell[[1L]]
        density <- cell[[2L]]
        fixture <- make_timing_fixture(width, density)
        if (mode == "verify") {
            verify_flip_cell(flip, other, fixture, width, density)
        } else {
            profile <- file.path(output, paste(width, density, "Rprofmem", sep = "."))
            measure_flip(mode, flip, fixture, profile, width, density)
        }
    }
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 11L)
sys.source(args[[2L]], envir = globalenv())
run_flip_benchmark(args)
