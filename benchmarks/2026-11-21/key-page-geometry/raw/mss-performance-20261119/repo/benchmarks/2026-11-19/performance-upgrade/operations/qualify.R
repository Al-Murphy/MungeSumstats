# Keep all test evidence before deciding whether the serial matrix can proceed.
qualify_suite <- function(tests, receipt, package = NULL) {
    if (file.exists(receipt) || !dir.create(receipt)) stop("Need a new qualification receipt directory")
    conditions <- list()
    record <- function(file, test, type, message) {
        conditions[[length(conditions) + 1L]] <<- data.frame(
            file = file, test = test, type = type, message = message, stringsAsFactors = FALSE)
    }
    write_tsv <- function(x, name) {
        utils::write.table(x, file.path(receipt, name), sep = "\t", quote = TRUE,
                           qmethod = "double", row.names = FALSE, fileEncoding = "UTF-8")
    }
    results <- tryCatch(withCallingHandlers(testthat::test_dir(
        tests, package = package, load_package = if (is.null(package)) "none" else "installed",
        reporter = "summary", stop_on_warning = FALSE, stop_on_failure = FALSE),
        warning = function(w) record("", "", "suite_warning", conditionMessage(w))),
        error = function(e) {
            record("", "", "suite_error", conditionMessage(e))
            NULL
        })
    counts <- data.frame(file=character(), test=character(), nb=integer(), passed=integer(),
                         skipped=logical(), failed=integer(), error=logical(), warning=integer())
    if (length(results)) {
        frame <- as.data.frame(results)
        counts <- frame[names(counts)]
        for (i in seq_along(results)) {
            for (event in results[[i]]$results) {
                if (!inherits(event, "expectation_success")) {
                    record(counts$file[i], counts$test[i], class(event)[1L], conditionMessage(event))
                }
            }
        }
    }
    write_tsv(counts, "tests.tsv")
    events <- if (length(conditions)) do.call(rbind, conditions) else
        data.frame(file=character(), test=character(), type=character(), message=character())
    write_tsv(events, "conditions.tsv")
    qualified <- nrow(counts) > 0L && sum(counts$passed) > 0L &&
        !any(counts$skipped | counts$failed > 0L | counts$error)
    write_tsv(data.frame(qualified=qualified, tests=nrow(counts), passed=sum(counts$passed),
                         skipped=sum(counts$skipped), failed=sum(counts$failed),
                         errors=sum(counts$error) + sum(events$type == "suite_error"),
                         test_warnings=sum(counts$warning), suite_warnings=sum(events$type == "suite_warning")),
              "qualification.tsv")
    if (!qualified) stop("Incomplete or failed testthat qualification; see ", receipt)
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    stopifnot(length(args) == 4L)
    .libPaths(c(args[[2]], args[[3]]), include.site = FALSE)
    qualify_suite(args[[1]], args[[4]], package = "MungeSumstats")
}
