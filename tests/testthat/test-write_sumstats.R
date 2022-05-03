test_that("write_sumstats works", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" ##&&
        ##.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        sumstats_dt <- MungeSumstats::formatted_example()
        
        run_tests <- function(sumstats_dt,
                              fileext,
                              tabix_index=FALSE,
                              write_vcf=FALSE,
                              standardise_headers=FALSE,
                              save_path_check=FALSE){
            message(fileext) 
            #### Check paths #### 
            message("=== write tests ===")
            path_out <- MungeSumstats::write_sumstats(
                sumstats_dt = sumstats_dt,
                save_path = tempfile(fileext = fileext),
                write_vcf = write_vcf,
                tabix_index = tabix_index,
                return_path = TRUE, 
                save_path_check = save_path_check
            )
            testthat::expect_true(file.exists(path_out))
            
            message("\n=== read tests ===")
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
                         fileext = ".tsv.bgz")
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
        # write_vcf=F
        ### Will write as tsv even though specified vcf suffix
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".vcf") 
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         fileext = ".vcf.gz") 
        
        #### VCF formats ####
        # write_vcf=T
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         write_vcf = TRUE,
                         fileext = ".vcf") 
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         write_vcf = TRUE,
                         fileext = ".vcf.gz") 
        ## with indexing
        # providing the correct suffix is important
       testthat::expect_failure(
           dat <- run_tests(sumstats_dt = sumstats_dt, 
                            write_vcf = TRUE,
                            tabix_index = TRUE,
                            fileext = ".vcf", 
                            save_path_check = FALSE) 
       )
       dat <- run_tests(sumstats_dt = sumstats_dt, 
                        write_vcf = TRUE,
                        tabix_index = TRUE,
                        fileext = ".vcf", 
                        save_path_check = TRUE) 
        dat <- run_tests(sumstats_dt = sumstats_dt, 
                         tabix_index = TRUE,
                         write_vcf = TRUE,
                         fileext = ".vcf.bgz") 
    }
})
