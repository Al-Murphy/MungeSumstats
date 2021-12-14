test_that("list_sumstats and parse_logs work", {
  
    #### Make some fake sums stats files ####
    # Let's pretend this is a different dataset each time.
    sumstats_path <- system.file(
        "extdata/ieu-a-298.tsv.gz",
        package = "MungeSumstats")
    #### Create fake log files ####
    # Let's pretend this is a different dataset each time.
    logs_path <- sumstats_path <- system.file(
        "extdata/MungeSumstats_log_msg.txt",
        package = "MungeSumstats")  
    #### iterate 10 files ####
    save_dir <- file.path(tempdir(),"data")
    n_files <- 10
    for(i in seq_len(n_files)){
        tmpdir <- file.path(save_dir,paste0("file",i))
        dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        tmpfile <- tempfile(tmpdir = tmpdir,
                            fileext = ".tsv.gz")
        print(tmpfile)
        #### sumstats ####
        file.copy(sumstats_path,
                  tmpfile)
        #### logs ####
        logs_path2 <- file.path(tmpdir,"logs")
        dir.create(logs_path2,showWarnings = FALSE)
        file.copy(logs_path,logs_path2)
    }
    #### List summary stats files ###
    munged_files <- MungeSumstats::list_sumstats(save_dir = save_dir)
    testthat::expect_length(munged_files,n_files)
    testthat::expect_length(unique(names(munged_files)), n_files)
    
    #### Parse log files ####
    log_data <- MungeSumstats::parse_logs(save_dir = save_dir)
    testthat::expect_true(methods::is(log_data,"data.frame"))
    testthat::expect_equal(nrow(log_data), n_files)
})
