test_that("compute_nsize works", {
  
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        file <- tempfile()
        sumstats_dt <- MungeSumstats::formatted_example(formatted = FALSE)
        ## pass mulitple N's - SNP level N not supported
        n_err <- 
          tryCatch(MungeSumstats::format_sumstats(sumstats_dt,
                                         ref_genome = "GRCh37",
                                         on_ref_genome = FALSE,
                                         strand_ambig_filter = FALSE,
                                         bi_allelic_filter = FALSE,
                                         allele_flip_check = FALSE,
                                         dbSNP=144,
                                         compute_n = c(5,5,7,8,10,12,12)
                                         ),
                   error = function(e) e,
                   warning = function(w) w
          )
        expect_true(is(n_err, "error"))
        
        #### Add N: numeric mode####
        sumstats_dt2 <- MungeSumstats::compute_nsize(sumstats_dt=sumstats_dt,
                                                     standardise_headers = TRUE, 
                                                     return_list = TRUE,
                                                     compute_n = 1000
                                                     )$sumstats_dt
        
        testthat::expect_false("N" %in% colnames(sumstats_dt))
        testthat::expect_true("N" %in% colnames(sumstats_dt2))
        
        #### Use existing N ####
        sumstats_dt3 <- compute_nsize(sumstats_dt = sumstats_dt2,
                                      compute_n = 500,
                                      return_list = FALSE)
        testthat::expect_equal(sumstats_dt3$N[1],1000)
        
        #### N_CAS/N_CON methods ####
        #### LDSC ####
        set.seed(1234)
        sumstats_dt$N_CAS <- as.integer(rnorm(nrow(sumstats_dt), mean = 500))
        sumstats_dt$N_CON <- as.integer(rnorm(nrow(sumstats_dt), mean = 1000))
        opts <- eval(as.list(args(compute_nsize))$compute_n)
        for(o in opts){
            print(paste0("Testing: compute_n='",o,"'"))
            sumstats_dt4 <- compute_nsize(sumstats_dt = sumstats_dt,
                                          compute_n = o,
                                          return_list = FALSE)
            if(o=="sum"){
                testthat::expect_true("N" %in% colnames(sumstats_dt4))
            } else {
                testthat::expect_true("Neff" %in% colnames(sumstats_dt4))
            } 
            remove(sumstats_dt4)
        } 
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})
