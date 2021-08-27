test_that("Check that imputation columns added correctly", {
    pth <- system.file("extdata", "eduAttainOkbay.txt",
        package = "MungeSumstats"
    )
    eduAttainOkbay <- data.table::fread(pth)
    # edit to make an rs id be imputed
    eduAttainOkbay[1, "MarkerName"] <-
        substring(
            eduAttainOkbay[1, "MarkerName"], 3,
            nchar(eduAttainOkbay[1, "MarkerName"])
        )
    # write to temp dir
    file <- tempfile()
    data.table::fwrite(eduAttainOkbay, file)

    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" &&
        .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        # run
        reformatted <- MungeSumstats::format_sumstats(file,
            ref_genome = "GRCh37",
            compute_z = TRUE,
            compute_n = 1001,
            ldsc_format = TRUE,
            imputation_ind = TRUE,
            allele_flip_check = TRUE
        )
        res <- data.table::fread(reformatted)
        col_headers <- names(res)
        imputat_cols <- c(
            col_headers[grepl("^IMPUTATION_", col_headers)],
            "flipped"["flipped" %in% col_headers],
            col_headers[grepl("^convert_", col_headers)]
        )
        # just check imputation columns exist
        expect_equal(length(imputat_cols) > 0, TRUE)
        # also check all have at least 1 value present
        have_value <- TRUE
        for (col_i in imputat_cols) {
            col_i_val <- res[[col_i]]
            if (length(col_i_val[!is.na(col_i_val)]) == 0) {
                  have_value <- FALSE
              }
        }
        expect_equal(have_value, TRUE)
    } else {
        expect_equal(is_32bit_windows, TRUE)
    }
})
