# Finite controls using real testthat and synthetic tests only. Retain each receipt.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !file.exists(args[[2]]), dir.create(args[[2]]))
source(file.path(args[[1]], "qualify.R"))
root <- args[[2]]
specs <- list(pass = "testthat::expect_identical(2L, 2L)",
              skip = 'testthat::skip("fixture dependency missing")',
              fail = 'testthat::expect_identical("actual", "expected")',
              warning = 'warning("fixture warning"); testthat::expect_true(TRUE)',
              error = 'stop("fixture error")')
for (kind in names(specs)) {
    suite <- file.path(root, kind)
    receipt <- file.path(root, paste0(kind, "-receipt"))
    stopifnot(dir.create(suite))
    writeLines(c("testthat::local_edition(3)",
                 'testthat::test_that("always passes", {testthat::expect_true(TRUE)})',
                 paste0('testthat::test_that("selected control", {', specs[[kind]], '})')),
               file.path(suite, "test-control.R"))
    passed_gate <- tryCatch(qualify_suite(suite, receipt), error=function(e) {
        if (!startsWith(conditionMessage(e), "Incomplete or failed testthat qualification")) stop(e)
        FALSE
    })
    summary <- read.delim(file.path(receipt, "qualification.tsv"))
    tests <- read.delim(file.path(receipt, "tests.tsv"))
    events <- read.delim(file.path(receipt, "conditions.tsv"))
    want <- kind %in% c("pass", "warning")
    stopifnot(identical(passed_gate, want), identical(summary$qualified, want),
              identical(tests$test, c("always passes", "selected control")),
              summary$skipped == as.integer(kind == "skip"),
              summary$failed == as.integer(kind == "fail"),
              summary$errors == as.integer(kind == "error"),
              summary$test_warnings == as.integer(kind == "warning"))
    if (kind != "pass") {
        stopifnot(any(events$file == "test-control.R" & events$test == "selected control" &
                      nzchar(events$message)))
    }
}
cat("Five real testthat qualification controls passed\n")
