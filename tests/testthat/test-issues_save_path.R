# Regression tests for GitHub issues on output paths (check_save_path) and reuse
# of saved files. Issues: #7, #10, #24, #33, #34, #49 (suffix part), #62.
run_check_save_path <- function(save_path, tabix_index = FALSE, ...) {
    MungeSumstats:::check_save_path(
        save_path = save_path,
        log_folder = tempdir(),
        log_folder_ind = FALSE,
        tabix_index = tabix_index,
        ...
    )
}

test_that("#7: tabix_index=TRUE rewrites a .gz save_path to .bgz", {
    f <- file.path(tempfile(), "ss.tsv.gz")
    expect_message(
        out <- run_check_save_path(f, tabix_index = TRUE),
        "Switching output to tabix-indexed format"
    )
    expect_equal(out$save_path, file.path(dirname(f), "ss.tsv.bgz"))
    expect_equal(out$extension, ".tsv.bgz")
})

test_that("#10: user-specified nested and relative save_paths are honoured", {
    sp <- file.path(tempfile(), "nested", "out.tsv.gz")
    out <- run_check_save_path(sp)
    expect_equal(out$save_path, sp)
    expect_true(dir.exists(dirname(sp)))
    expect_equal(out$sep, "\t")
    expect_equal(out$file_type, ".tsv.gz")
    # "./" is expanded against the working directory
    wd <- tempfile("wd")
    dir.create(wd)
    old <- setwd(wd)
    on.exit(setwd(old), add = TRUE)
    wd_now <- getwd()
    out_rel <- run_check_save_path("./rel_out/out.tsv.gz")
    setwd(old)
    expect_equal(out_rel$save_path, paste0(wd_now, "/rel_out/out.tsv.gz"))
    expect_true(dir.exists(file.path(wd, "rel_out")))
})

test_that("#33: saving to tempdir() prints NOTEs; other folders do not", {
    msgs <- testthat::capture_messages(
        MungeSumstats:::check_save_path(
            save_path = file.path(tempdir(), "a.tsv.gz"),
            log_folder = tempdir(),
            log_folder_ind = TRUE,
            tabix_index = FALSE
        )
    )
    msgs <- paste(msgs, collapse = "")
    expect_match(msgs, "Formatted results will be saved to `tempdir\\(\\)`")
    expect_match(msgs, "Log results will be saved to `tempdir\\(\\)`")
    msgs2 <- testthat::capture_messages(
        run_check_save_path(file.path(tempfile(), "a.tsv.gz"))
    )
    expect_false(any(grepl("tempdir()", msgs2, fixed = TRUE)))
})

test_that("#33: '~' in save_path is expanded to the home directory", {
    home <- tempfile("home")
    dir.create(home)
    old_home <- Sys.getenv("HOME")
    Sys.setenv(HOME = home)
    on.exit(Sys.setenv(HOME = old_home), add = TRUE)
    # never create directories in the real HOME
    skip_if_not(path.expand("~") == home, "'~' is not expanded from HOME")
    out <- run_check_save_path("~/x_mss_test/a.tsv.gz")
    Sys.setenv(HOME = old_home)
    expect_false(startsWith(out$save_path, "~"))
    expect_equal(out$save_path, file.path(home, "x_mss_test", "a.tsv.gz"))
    expect_true(dir.exists(file.path(home, "x_mss_test")))
})

test_that("#34: the save_path suffix follows write_vcf", {
    td <- tempfile()
    # VCF suffix but write_vcf=FALSE -> tabular .tsv.gz
    for (f in c("x.vcf", "x.vcf.gz", "x.vcf.bgz")) {
        expect_message(
            out <- run_check_save_path(file.path(td, f), write_vcf = FALSE),
            "save_path suggests VCF output but write_vcf=FALSE"
        )
        expect_equal(out$save_path, file.path(td, "x.tsv.gz"))
        expect_equal(out$file_type, ".tsv")
        expect_equal(out$extension, ".tsv.gz")
    }
    # tabular suffix but write_vcf=TRUE -> .vcf.bgz
    out <- run_check_save_path(file.path(td, "x.tsv.gz"), write_vcf = TRUE)
    expect_equal(out$save_path, file.path(td, "x.vcf.bgz"))
    expect_equal(out$file_type, "vcf")
})

test_that("#34: only the file suffix is rewritten, not directory names", {
    td <- tempfile()
    out <- run_check_save_path(file.path(td, "my.txt.dir", "out.tsv.gz"),
                               write_vcf = TRUE)
    expect_equal(out$save_path, file.path(td, "my.txt.dir", "out.vcf.bgz"))
    out <- run_check_save_path(file.path(td, "my.vcf.dir", "out.vcf.gz"),
                               write_vcf = FALSE)
    expect_equal(out$save_path, file.path(td, "my.vcf.dir", "out.tsv.gz"))
})

test_that("#49: output format is matched on the end of save_path", {
    td <- tempfile()
    for (ext in c(".tsv", ".tsv.gz", ".csv.gz", ".txt.bgz")) {
        out <- run_check_save_path(file.path(td, paste0("ss", ext)),
                                   verbose = FALSE)
        expect_equal(out$extension, ext)
        expect_equal(out$sep, if (startsWith(ext, ".csv")) "," else "\t")
    }
    expect_error(
        run_check_save_path(file.path(td, "ss.parquet"), verbose = FALSE),
        "save_path file format not recognized"
    )
})

test_that("#62: a save_path folder that cannot be created fails early", {
    # the parent is a regular file, so no OS can create the folder
    f <- tempfile()
    file.create(f)
    expect_error(
        run_check_save_path(file.path(f, "sub", "out.tsv.gz"), verbose = FALSE),
        "Unable to create folder for save_path"
    )
})

test_that("#24: a previously formatted save_path short-circuits format_sumstats", {
    dt <- data.table::fread(
        system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats"),
        nThread = 1
    )
    sp <- file.path(tempfile(), "prev.tsv.gz")
    dir.create(dirname(sp))
    MungeSumstats::write_sumstats(sumstats_dt = dt, save_path = sp)
    # a path that does not exist proves the munging pipeline was skipped
    expect_message(
        out <- MungeSumstats::format_sumstats(path = "does/not/exist.tsv",
                                              save_path = sp,
                                              force_new = FALSE),
        "Importing previously formatted file"
    )
    expect_equal(out, sp)
    res <- MungeSumstats::format_sumstats(path = "does/not/exist.tsv",
                                          save_path = sp,
                                          return_data = TRUE)
    expect_true(data.table::is.data.table(res))
    expect_equal(dim(res), c(93, 9))
    expect_equal(names(res), names(dt))
})

test_that("#24: a global `sumstats_return` does not replace the cached data", {
    skip_if(exists("sumstats_return", envir = globalenv(), inherits = FALSE))
    dt <- data.table::fread(
        system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats"),
        nThread = 1
    )
    sp <- file.path(tempfile(), "prev.tsv.gz")
    dir.create(dirname(sp))
    MungeSumstats::write_sumstats(sumstats_dt = dt, save_path = sp)
    assign("sumstats_return", list(sumstats_dt = "WRONG"), envir = globalenv())
    on.exit(rm("sumstats_return", envir = globalenv()), add = TRUE)
    res <- MungeSumstats::format_sumstats(path = "does/not/exist.tsv",
                                          save_path = sp,
                                          return_data = TRUE)
    expect_true(data.table::is.data.table(res))
    expect_equal(nrow(res), 93)
})
