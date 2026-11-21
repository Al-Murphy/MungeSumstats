# Regression tests for GitHub issues on run logs (format_sumstats log files and sinks,
# report_summary, parse_logs). Issues: #4/#32, #61/#131/#207/#80, #68.
fast_ref_args <- list(
    ref_genome = "GRCh37", dbSNP = 144, on_ref_genome = FALSE,
    strand_ambig_filter = FALSE, bi_allelic_filter = FALSE,
    allele_flip_check = FALSE
)
skip_if_no_grch37 <- function() {
    # loading the SNPlocs namespace warns about replaced imports
    suppressWarnings({
        skip_if_not_installed("SNPlocs.Hsapiens.dbSNP144.GRCh37")
        skip_if_not_installed("BSgenome.Hsapiens.1000genomes.hs37d5")
    })
}
log_is_open <- function(name) {
    any(grepl(name, showConnections(all = TRUE)[, "description"],
              fixed = TRUE))
}

test_that("#68: parse_logs finds current '<id>_log_msg.txt' and old log names", {
    root <- tempfile("parse_logs")
    for (id in c("newid", "oldid")) {
        dir.create(file.path(root, id, "logs"), recursive = TRUE)
    }
    ex <- MungeSumstats:::logs_example()
    file.copy(ex, file.path(root, "newid", "logs", "newid_log_msg.txt"))
    file.copy(ex, file.path(root, "oldid", "logs", "MungeSumstats_log_msg.txt"))
    # output logs are not message logs
    file.copy(ex, file.path(root, "newid", "logs", "newid_log_output.txt"))
    res <- MungeSumstats::parse_logs(save_dir = root)
    expect_equal(nrow(res), 2)
    expect_setequal(res$id, c("newid", "oldid"))
    new <- res[res$id == "newid", ]
    expect_equal(basename(new$log_path), "newid_log_msg.txt")
    expect_equal(new$rows_start, 93)
    expect_equal(new$snps_start, 93)
    expect_equal(new$sig_snps_start, 64)
    expect_equal(new$chroms_start, 20)
    expect_equal(new$rows_end, 93)
})

test_that("#4/#32: report_summary wording round-trips through parse_report", {
    dt <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs2", "rs3"),
        CHR = c("1", "1", "2", "X"),
        P = c(1e-9, 0.5, 1e-10, 0.2)
    )
    msg <- testthat::capture_messages(
        MungeSumstats:::report_summary(sumstats_dt = data.table::copy(dt),
                                       orig_dims = c(8, 3))
    )
    l <- strsplit(paste(msg, collapse = ""), "\n")[[1]]
    expect_equal(l, c(
        "Summary statistics report:",
        "   - 4 rows (50% of original 8 rows)",
        "   - 3 unique variants",
        "   - 2 genome-wide significant variants (P<5e-8)",
        "   - 3 chromosomes"
    ))
    parsed <- vapply(1:4, function(i) {
        MungeSumstats:::parse_report(l = l, entry = 1, line = i)
    }, integer(1))
    expect_equal(parsed, c(4L, 3L, 2L, 3L))
    # thousands separators are stripped when parsed
    big <- testthat::capture_messages(
        MungeSumstats:::report_summary(
            sumstats_dt = data.table::data.table(BP = seq_len(1500))
        )
    )
    l_big <- strsplit(paste(big, collapse = ""), "\n")[[1]]
    expect_equal(l_big[2], "   - 1,500 rows")
    expect_equal(MungeSumstats:::parse_report(l = l_big, entry = 1, line = 1),
                 1500L)
})

test_that("#61/#131/#207/#80: message logs are named from save_path and sinks are restored", {
    skip_if_no_grch37()
    edu <- system.file("extdata", "eduAttainOkbay.txt",
                       package = "MungeSumstats")
    td <- tempfile("logs131")
    dir.create(td)
    ld <- file.path(td, "logs")
    n0 <- sink.number()
    n_msg <- sink.number(type = "message")
    on.exit(while (sink.number() > n0) sink(), add = TRUE)
    #### success path, inside a caller's own output sink (knitr-like) ####
    sink(tempfile())
    n_caller <- sink.number()
    warns <- character(0)
    res <- withCallingHandlers(
        do.call(MungeSumstats::format_sumstats, c(fast_ref_args, list(
            path = edu,
            save_path = file.path(td, "gwasA.tsv.gz"),
            log_folder = ld,
            log_mungesumstats_msgs = TRUE,
            force_new = TRUE
        ))),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    n_after <- sink.number()
    if (sink.number() > n0) sink()
    # #207/#80: the caller's sink survives, no double unsink, log closed
    expect_equal(n_after, n_caller)
    expect_equal(sink.number(type = "message"), n_msg)
    expect_false(any(grepl("no sink to remove", warns)))
    expect_false(log_is_open("gwasA_log_msg.txt"))
    # #131: log files are named after the save_path stem
    expect_equal(basename(res$log_files$MungeSumstats_log_msg),
                 "gwasA_log_msg.txt")
    expect_equal(basename(res$log_files$MungeSumstats_log_output),
                 "gwasA_log_output.txt")
    expect_true(all(file.exists(unlist(res$log_files))))
    # #61: messages are recorded in the log
    log_a <- readLines(res$log_files$MungeSumstats_log_msg)
    expect_true(any(grepl("Summary statistics report:", log_a, fixed = TRUE)))
    expect_true(any(grepl("Done munging in", log_a, fixed = TRUE)))
    #### error path, with a second save_path ####
    bad <- tempfile(fileext = ".foo")
    writeLines("garbage", bad)
    expect_error(
        do.call(MungeSumstats::format_sumstats, c(fast_ref_args, list(
            path = bad,
            save_path = file.path(td, "gwasB.tsv.gz"),
            log_folder = ld,
            log_mungesumstats_msgs = TRUE
        ))),
        "Unrecognized file format"
    )
    expect_equal(sink.number(), n0)
    expect_equal(sink.number(type = "message"), n_msg)
    expect_false(log_is_open("gwasB_log_msg.txt"))
    # #131: the second run gets its own logs and leaves the first untouched
    expect_true(file.exists(file.path(ld, "gwasB_log_msg.txt")))
    expect_equal(readLines(file.path(ld, "gwasA_log_msg.txt")), log_a)
})

test_that("#80: rerunning with force_new=TRUE does not append to old logs", {
    skip_if_no_grch37()
    edu <- system.file("extdata", "eduAttainOkbay.txt",
                       package = "MungeSumstats")
    td <- tempfile("logs80")
    dir.create(td)
    for (i in 1:2) {
        res <- do.call(MungeSumstats::format_sumstats, c(fast_ref_args, list(
            path = edu,
            save_path = file.path(td, "run.tsv.gz"),
            log_folder = file.path(td, "logs"),
            log_mungesumstats_msgs = TRUE,
            force_new = TRUE
        )))
    }
    log <- readLines(res$log_files$MungeSumstats_log_msg)
    expect_equal(sum(grepl("Done munging in", log, fixed = TRUE)), 1)
})

test_that("#207: a caller's message sink is restored after format_sumstats", {
    skip_if_no_grch37()
    td <- tempfile("logs207")
    dir.create(td)
    bad <- tempfile(fileext = ".foo")
    writeLines("garbage", bad)
    con <- file(tempfile(), open = "wt")
    sink(con, type = "message")
    n_msg <- sink.number(type = "message")
    try(do.call(MungeSumstats::format_sumstats, c(fast_ref_args, list(
        path = bad,
        save_path = file.path(td, "x.tsv.gz"),
        log_folder = file.path(td, "logs"),
        log_mungesumstats_msgs = TRUE
    ))), silent = TRUE)
    n_after <- sink.number(type = "message")
    sink(type = "message")
    close(con)
    expect_equal(n_after, n_msg)
})
