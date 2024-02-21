test_that("check_info_score works", {
    
    sumstats_dt <- MungeSumstats::formatted_example()  
    ### Make fake INFO score ####
    sumstats_dt$INFO <- rep(c(1,.95,.8,0), ceiling(nrow(sumstats_dt)/4))[
        seq_len(nrow(sumstats_dt))
    ]
    #### Make fake check_save_out list ####
    check_save_out <- list(sep="\t",
                           extension=".tsv.gz",
                           log_folder=tempdir())
    log_files <- list(
        MungeSumstats_log_msg=file.path(tempdir(),
                                        "MungeSumstats_log_msg.txt"),
        MungeSumstats_log_output=file.path(tempdir(),
                                           "MungeSumstats_log_output.txt"))
    
    #### INFO_filter=.9  ####
    sumstats <- MungeSumstats:::check_info_score(sumstats_dt=sumstats_dt, 
                                                 INFO_filter=.9,
                                                 log_folder_ind=TRUE,
                                                 check_save_out=check_save_out, 
                                                 tabix_index=FALSE, 
                                                 nThread=1, 
                                                 log_files=log_files) 
    testthat::expect_equal(names(sumstats),c( "sumstats_dt","log_files" )) 
    testthat::expect_equal(nrow(sumstats$sumstats_dt),47)
    testthat::expect_lt(nrow(sumstats$sumstats_dt), nrow(sumstats_dt))
    dropped_snps <- data.table::fread(sumstats$log_files$info_filter)
    testthat::expect_equal(nrow(sumstats_dt) - nrow(dropped_snps), 
                           nrow(sumstats$sumstats_dt))
    
    
    #### INFO_filter=0  ####
    sumstats <- MungeSumstats:::check_info_score(sumstats_dt=sumstats_dt, 
                                                 INFO_filter=0,
                                                 log_folder_ind=TRUE,
                                                 check_save_out=check_save_out, 
                                                 tabix_index=FALSE, 
                                                 nThread=1, 
                                                 log_files=log_files) 
    testthat::expect_equal(names(sumstats),c( "sumstats_dt","log_files" )) 
    testthat::expect_equal(sumstats$sumstats_dt,sumstats_dt)  
    
    #### All INFO==1  ####
    sumstats_dt$INFO <- 1
    sumstats <- MungeSumstats::: check_info_score(sumstats_dt=sumstats_dt, 
                                                  INFO_filter=.9,
                                                  log_folder_ind=TRUE,
                                                  check_save_out=check_save_out, 
                                                  tabix_index=FALSE, 
                                                  nThread=1, 
                                                  log_files=log_files) 
    testthat::expect_equal(names(sumstats),c( "sumstats_dt","log_files" )) 
    testthat::expect_equal(sumstats$sumstats_dt,sumstats_dt)  
})