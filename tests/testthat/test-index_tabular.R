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
    }    
    else{
      expect_equal(is_32bit_windows, TRUE)
    }
})
