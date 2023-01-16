test_that("index_tabular works", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
  
    #### Don't run on Windows at all for now, since tabix-
    is_32bit_windows <-
      .Platform$OS.type == "windows"## && .Platform$r_arch == "i386"
    if (!is_32bit_windows) { 
        
      sumstats_dt <- MungeSumstats::formatted_example() 
      path <- tempfile(fileext = ".tsv")
      MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt,
                                    save_path = path)
      #### Index file ####
      tbx_file <- MungeSumstats::index_tabular(path = path) 
      
      #### Test that file exists ####
      testthat::expect_true(endsWith(tbx_file,".bgz")) 
      testthat::expect_true(file.exists(tbx_file)) 
      
      #### Test that header exists ###
      header <- MungeSumstats::read_header(path = tbx_file)
      testthat::expect_equal(ncol(header), 9)
      #### Test that all chr1 SNPs are present ####
      tabixRange <- paste0(1,":",
                           min(sumstats_dt$BP),"-",
                           max(sumstats_dt$BP))
      query <-  MungeSumstats::read_header(path = tbx_file, 
                                           n = NULL)
      query[,CHR:=as.character(CHR)] 
      testthat::expect_equal(query,sumstats_dt) 
      
      #### --- Test sorting --- ####
      #### GenomicRanges ####
      query_unsorted <- data.table::copy(query)
      data.table::setorderv(x = query_unsorted, 
                            cols = c("BETA"))
      query_resorted <- MungeSumstats:::sort_coords(
          sumstats_dt = data.table::copy(query_unsorted), 
          sort_method = "GenomicRanges")
      testthat::expect_equal(all.equal(query,query_resorted,
                                       ignore.col.order = TRUE),TRUE)
      testthat::expect_failure(
          testthat::expect_equal(query$SNP,query_unsorted$SNP)
      )
      #### data.table ####
      query_unsorted <- data.table::copy(query)
      data.table::setorderv(x = query_unsorted, 
                            cols = c("BETA"))
      query_resorted <- MungeSumstats:::sort_coords(
          sumstats_dt = data.table::copy(query_unsorted), 
          sort_method = "data.table")
      testthat::expect_equal(query,query_resorted)
      testthat::expect_failure(
          testthat::expect_equal(query$SNP,query_unsorted$SNP)
      )
    }    
    else{
      expect_equal(is_32bit_windows, TRUE)
    }
})
