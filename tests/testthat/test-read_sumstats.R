test_that("Test that different kinds of files can be read in.", {
    ## Call uses reference genome as default with more than 2GB of memory,
    ## which is more than what 32-bit Windows can handle so remove tests
    is_32bit_windows <-
        .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        ### read data
        path <- system.file("extdata", "eduAttainOkbay.txt", 
                            package = "MungeSumstats")
        sumstats_dt <- data.table::fread(path, nThread = 1)
        sumstats_dt <- standardise_header(
            sumstats_dt = sumstats_dt,
            return_list = FALSE,
        )
        col_order <- colnames(sumstats_dt)
        ### Test tsv
        tsv_ss <- sumstats_dt
        ### Test csv
        path_csv <- tempfile()
        data.table::fwrite(sumstats_dt, path_csv, sep = ",")
        csv_ss <- MungeSumstats::read_sumstats(path = path_csv, 
                                               standardise_headers = TRUE)
        #### Test space-separated
        path_space <- tempfile()
        data.table::fwrite(sumstats_dt, path_space, sep = " ")
        space_ss <- MungeSumstats::read_sumstats(path = path_space,
                                                 standardise_headers = TRUE)
        ### Test VCF
        vcf_tmp <- tempfile(fileext = ".vcf.gz")
        vcf_tmp <- write_sumstats(
            sumstats_dt = sumstats_dt,
            save_path = vcf_tmp,
            write_vcf = TRUE,
            return_path = TRUE,
            save_path_check = TRUE
        )
        vcf_ss <- MungeSumstats::read_sumstats(path = vcf_tmp,
                                               standardise_headers = TRUE) 
        for (x in c("FRQ", "BETA", "SE", "P")) {
            vcf_ss[, (x) := as.numeric(get(x))]
        }
        vcf_ss[,c("ID","END"):=NULL]
        data.table::setcolorder(vcf_ss, col_order) 
        
        #### Test that all dts are the same
        objs <- mget(c("tsv_ss", "csv_ss", "space_ss", "vcf_ss"))
        #### Ensure all keys are the same ####
        for(obj in objs){
            data.table::setkey(obj,SNP)
        }
        iterative_tests <- outer(objs, objs, Vectorize(all.equal))
        if (any(is(unlist(iterative_tests[1, ]), "character"))) {
            all_same <- FALSE
        } else {
            all_same <- all(colSums(iterative_tests) == nrow(iterative_tests))
        } 
        testthat::expect_equal(all_same, TRUE)
    }    
    else{
        expect_equal(is_32bit_windows, TRUE)
    }
})