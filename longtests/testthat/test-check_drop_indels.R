test_that("check_range_p_val works", {
    
    sumstats_dt <- MungeSumstats:::formatted_example()
    #copy the org
    sumstats_dt_org <- copy(sumstats_dt)
    #add Indel
    sumstats_dt_indel <- sumstats_dt[1,]
    sumstats_dt_indel[1,A2:='CG']
    sumstats_dt <- data.table::rbindlist(list(sumstats_dt,
                                              sumstats_dt_indel))
    check_save_out <- check_save_path(
      save_path = tempfile(fileext = ".tsv.gz"),
      log_folder = tempdir(),
      log_folder_ind = TRUE,
      tabix_index = FALSE,
      write_vcf = FALSE
    )
    
    sumstats <- check_drop_indels(sumstats_dt = sumstats_dt,
                                  drop_indels = TRUE,
                                  path = tempfile(fileext = ".tsv.gz"),
                                  log_folder_ind = TRUE,
                                  check_save_out=check_save_out,
                                  tabix_index = FALSE,
                                  nThread=1,
                                  log_files=vector(mode = "list"))
    testthat::expect_equal(
        sumstats$sumstats_dt,sumstats_dt_org)
})
