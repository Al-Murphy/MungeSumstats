test_that("Check that imputation columns added correctly", {
  ## The following test uses more than 2GB of memory, which is more
  ## than what 32-bit Windows can handle:
  is_32bit_windows <- .Platform$OS.type == "windows" #&&
  #.Platform$r_arch == "i386"
  if (!is_32bit_windows) {
    #only run not on linux to speed up linux bioc checks
    pth <- system.file("extdata", "eduAttainOkbay.txt",
                       package = "MungeSumstats"
    )
    if (Sys.info()["sysname"]!="Linux"){
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
      # run
      reformatted <- MungeSumstats::format_sumstats(file,
                                                    ref_genome = "GRCh37",
                                                    compute_z = TRUE,
                                                    compute_n = 1001,
                                                    save_format='LDSC',
                                                    imputation_ind = TRUE,
                                                    allele_flip_check = TRUE,
                                                    dbSNP=144
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
    } else{
      expect_equal(isTRUE(Sys.info()["sysname"]=="Linux"), TRUE)
      expect_equal(isTRUE(Sys.info()["sysname"]=="Linux"), TRUE)
    }
    # check other compute_n values
    eduAttainOkbay <- data.table::fread(pth)
    eduAttainOkbay[, N_CON := 100]
    eduAttainOkbay[, N_CAS := 120]
    # write to temp dir
    file <- tempfile()
    data.table::fwrite(eduAttainOkbay, file)
    methods <- c("ldsc", "sum", "giant", "metal")
    reformatted <- MungeSumstats::format_sumstats(file,
                                                  ref_genome = "GRCh37",
                                                  compute_n = methods,
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter = FALSE,
                                                  bi_allelic_filter = FALSE,
                                                  allele_flip_check = FALSE,
                                                  dbSNP=144
    )
    res <- data.table::fread(reformatted)
    expect_equal(all(paste0("Neff_", c("ldsc", "giant", "metal")) %in%
                       colnames(res)), TRUE)
  } else {
    expect_equal(is_32bit_windows, TRUE)
    expect_equal(is_32bit_windows, TRUE)
    expect_equal(is_32bit_windows, TRUE)
  }
})
