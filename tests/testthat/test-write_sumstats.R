test_that("write_sumstats works", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
        ##.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        sumstats_dt <- MungeSumstats:::formatted_example()
        
        
        run_tests <- function(sumstats_dt,
                              fileext,
                              tabix_index=FALSE,
                              write_vcf=FALSE,
                              standardise_headers=FALSE){
            message(fileext) 
            #### Check paths ####
            path_in <- tempfile(fileext = fileext)
            check <- MungeSumstats:::check_save_path(save_path = path_in, 
                                                     log_folder = tempdir(),
                                                     log_folder_ind = FALSE,
                                                     tabix_index = tabix_index, 
                                                     write_vcf = write_vcf)
            path_in <- check$save_path
            
            path_out <- MungeSumstats:: write_sumstats(
                sumstats_dt = sumstats_dt,
                save_path = path_in,
                write_vcf = write_vcf,
                tabix_index = tabix_index,
                return_path = TRUE
            )
            testthat::expect_true(file.exists(path_out))
            dat <- MungeSumstats::read_sumstats(
                path = path_out, 
                standardise_headers = standardise_headers)
            testthat::expect_equal(nrow(dat), nrow(sumstats_dt)) 
            return(dat)
        }
        
        #### Tabular formats ####
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".tsv")
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".tsv.gz")
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".tsv.gz", 
                         tabix_index = TRUE)
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".tsv.bgz", 
                         tabix_index = TRUE)
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".csv")
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".csv.gz")
        
        #### VCF formats ####
        # write_vcf=F
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".vcf") 
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".vcf.gz") 
        # write_vcf=T
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         write_vcf = TRUE,
                         fileext = ".vcf") 
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         write_vcf = TRUE,
                         tabix_index = TRUE,
                         fileext = ".vcf") 
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         tabix_index = TRUE,
                         write_vcf = TRUE,
                         fileext = ".vcf.bgz") 
    }
})
