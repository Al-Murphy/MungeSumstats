test_that("Check that a file is in tabular format", {
    ### read data
    path <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
    sumstats_dt =  data.table::fread(path, nThread = 1)
    ### Test tsv 
    header <- read_header(path = path) 
    is_tabular_tsv <- check_tabular(header = header)
    ### Test csv 
    path_csv <- tempfile()
    data.table::fwrite(sumstats_dt, path_csv, sep = ",")
    header <- read_header(path = path_csv) 
    is_tabular_csv <- check_tabular(header = header) 
    #### Test space-separated  
    path_space <- tempfile()
    data.table::fwrite(sumstats_dt, path_space, sep = " ")
    header <- read_header(path = path_space) 
    is_tabular_space <- check_tabular(header = header)   
  expect_equal(all.equal(is_tabular_tsv, is_tabular_csv, is_tabular_space), TRUE)
})
