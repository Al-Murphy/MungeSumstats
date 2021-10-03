test_that("index_tabular works", {
  
    eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
                                     package = "MungeSumstats")
    sumstats_dt <- data.table::fread(eduAttainOkbayPth, nThread = 1)
    sumstats_dt <-
    MungeSumstats:::standardise_sumstats_column_headers_crossplatform(
        sumstats_dt = sumstats_dt)$sumstats_dt
    sumstats_dt <- MungeSumstats:::sort_coords(sumstats_dt = sumstats_dt)
    path <- tempfile(fileext = ".tsv")
    MungeSumstats::write_sumstats(sumstats_dt = sumstats_dt,
                                  save_path = path)
    #### Index file ####
    indexed_file <- MungeSumstats::index_tabular(path = path) 
    
    #### Test that file exists ####
    testthat::expect_true(endsWith(indexed_file,".bgz"))
    testthat::expect_true(file.exists(indexed_file))
    
    #### Test that header exists ###
    header <- seqminer::tabix.read.header(tabixFile = indexed_file)$header
    header <- colnames(data.table::fread(text = rep(header,2)))
    testthat::expect_length(header, 9)
    #### Test that all chr1 SNPs are present ####
    tabixRange <- paste0(1,":",
                         min(sumstats_dt$BP),"-",
                         max(sumstats_dt$BP))
    query <- seqminer::tabix.read.table(tabixFile = indexed_file,
                                        tabixRange = tabixRange) 
    # query <- data.table::data.table(query)
    testthat::expect_equal(query$SNP, subset(sumstats_dt,CHR==1)$SNP) 
})
