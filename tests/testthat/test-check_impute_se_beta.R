test_that("Check that imputation columns added correctly", {
    ## The following test uses more than 2GB of memory, which is more
    ## than what 32-bit Windows can handle:
    is_32bit_windows <- .Platform$OS.type == "windows" #&&
    #.Platform$r_arch == "i386"
    if (!is_32bit_windows) {
        #Our VCF file contains a BETA and SE column, let's impute these and 
        #compare the difference
        # write the ALS GWAS, VCF file to a temp file for testing
        vcf_path <- system.file("extdata", "ALSvcf.vcf", 
                                package = "MungeSumstats")
        sumstats_dt <- MungeSumstats:::read_vcf(path = vcf_path)
        #add in known sample size of study
        sumstats_dt[,N:=293723]
        # write to temp dir
        file <- tempfile()
        data.table::fwrite(sumstats_dt, file)
        #run control through
        control <- MungeSumstats::format_sumstats(file,
                                                  ref_genome = "GRCh37",
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter = FALSE,
                                                  bi_allelic_filter = FALSE,
                                                  allele_flip_check = FALSE,
                                                  imputation_ind = TRUE,
                                                  log_folder_ind = TRUE,
                                                  INFO_filter = 0,
                                                  dbSNP=144
        )
        control_dt <- data.table::fread(control$sumstats)
        #now impute SE
        sumstats_dt_org <- copy(sumstats_dt)
        sumstats_dt[,SE:=NULL]
        data.table::fwrite(sumstats_dt, file)
        reformatted <- MungeSumstats::format_sumstats(file,
                                                      ref_genome = "GRCh37",
                                                      impute_se = TRUE,
                                                      on_ref_genome = FALSE,
                                                      strand_ambig_filter =
                                                          FALSE,
                                                      bi_allelic_filter = FALSE,
                                                      allele_flip_check = FALSE,
                                                      imputation_ind = TRUE,
                                                      log_folder_ind = TRUE,
                                                      INFO_filter = 0,
                                                      dbSNP=144
        )
        #first check it was imputed
        imp_se_dt <- data.table::fread(reformatted$sumstats)
        expect_equal("SE" %in% colnames(imp_se_dt), TRUE)
        expect_equal("IMPUTATION_SE" %in% colnames(imp_se_dt), TRUE)
        #check how close it is to actual
        setkey(imp_se_dt,SNP)
        setkey(control_dt,SNP)
        control_dt[imp_se_dt,imp_se:=i.SE]
        #check diff as percentage of value
        control_dt[,diff_se:=abs((SE-imp_se)/SE)]
        #expect a mean difference below 0.1%
        testthat::expect_equal(mean(control_dt$diff_se)<0.001, TRUE)
        
        
        #-----------
        # Now let's impute beta
        # use sumstats with Z and beta:
        # Genome-wide association meta-analysis (N=269,867) identifies new 
        # genetic and functional links to intelligence.
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6411041/
        # PMID: 29942086 PMCID: PMC6411041 DOI: 10.1038/s41588-018-0152-6
        #manually writing 25 rows because of package size limits
        sumstats_dt_org <-
            data.table::data.table(SNP=c('rs12184267','rs12184277','rs12184279',
                                         'rs116801199','rs12565286','rs2977670',
                                         'rs28454925','rs116720794','rs4951859',
                                         'rs146277091','rs3094315','rs3131972',
                                         'rs3115860','rs2073813','rs12184312',
                                         'rs12184325','rs3131969','rs3131968',
                                         'rs12184313','rs3131967','rs3115859',
                                         'rs10454459','rs3131966','rs3131965',
                                         'rs3115858'),
                                   CHR=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                         1,1,1,1,1,1),
                                   BP=c(715265,715367,717485,720381,721290,
                                        723891,
                                        726794,729632,729679,752478,752566,
                                        752721,
                                        753405,753541,754063,754105,754182,
                                        754192,
                                        754211,754334,754503,754629,754964,
                                        755775,
                                        755890),
                                   A1=c('C','A','C','G','G','G','C','C','C',
                                        'G','G',
                                        'A','C','G','G','C','A','A','G','T',
                                        'G','A','C','A','A'),
                                   A2=c('T','G','A','T','C','C','G','T','G','A',
                                        'A',
                                        'G','A','A','T','T','G','G','A','C','A',
                                        'G','T','G','T'),
                                   FRQ=c(0.9591931,0.9589313,0.9594241,0.957838,
                                         0.9576224,0.06312,0.9590545,0.9589005,
                                         0.183523,0.9592701,0.180135,0.186495,
                                         0.147228,0.848645,0.9590545,0.9590083,
                                         0.154419,0.154389,0.9587773,0.153927,
                                         0.185618,0.9585925,0.184909,0.177333,
                                         0.147844),
                                   Z=c(-0.916,-0.656,-1.05,-0.3,-0.566,0.253,
                                       -0.539,
                                       -0.27,0.208,-0.263,0.351,0.61,0.48,
                                       -0.441,
                                       -0.181,-0.217,0.344,0.434,-0.212,0.509,
                                       0.823,-0.351,0.828,0.808,0.414),
                                   BETA=c(-0.0068872979,-0.0049144905,
                                          -0.0079116035,
                                          -0.0022174032,-0.0041742154,
                                          0.0015498434,
                                          -0.0040387267,-0.0020185581,
                                          0.0008058415,
                                          -0.0019709245,0.0013480759,
                                          0.0023117892,
                                          0.0019985366,-0.0018164389,
                                          -0.0013546785,
                                          -0.0016223525,0.0013769675,
                                          0.001760705,
                                          -0.0015823145,0.0020819201,
                                          0.003133556,
                                          -0.0026146926,0.003159564,
                                          0.0032010789,0.0017171769),
                                   SE=c(0.007518884,0.007491601,0.00753486,
                                        0.007391344,0.007374939,0.006125863,
                                        0.007492999,0.007476141,0.003874238,
                                        0.007494009,0.003840672,0.003789818,
                                        0.004163618,0.004118909,0.007484411,
                                        0.007476279,0.004002812,0.004056924,
                                        0.007463748,0.004090216,0.00380748,
                                        0.007449267,0.003815899,0.003961731,
                                        0.00414777),
                                   P=c(0.3598,0.5116,0.2938,0.7644,0.5711,
                                       0.8006,
                                       0.5896,0.7872,0.835,0.7924,0.7257,0.5421,
                                       0.631,0.6592,0.8567,0.8281,0.7309,0.6646,
                                       0.8322,0.6105,0.4103,0.7257,0.4074,0.419,
                                       0.6787),
                                   N=c(225955,226215,226224,226626,226528,
                                       225312,226782,226989,222312,227870,
                                       229517,229459,229723,229447,227303,
                                       227552,238992,232696,227092,229485,
                                       228163,227002,227830,218366,230684))
        
        sumstats_dt <- copy(sumstats_dt_org)
        sumstats_dt[,BETA:=NULL]
        data.table::fwrite(sumstats_dt, file)
        reformatted <- MungeSumstats::format_sumstats(file,
                                                      ref_genome = "GRCh37",
                                                      impute_beta = TRUE,
                                                      on_ref_genome = FALSE,
                                                      strand_ambig_filter = 
                                                          FALSE,
                                                      bi_allelic_filter = FALSE,
                                                      allele_flip_check = FALSE,
                                                      imputation_ind = TRUE,
                                                      log_folder_ind = TRUE,
                                                      dbSNP=144
        )
        imp_beta_dt <- data.table::fread(reformatted$sumstats)
        #first check it was imputed
        testthat::expect_equal("BETA" %in% colnames(imp_beta_dt), TRUE)
        #check imputation column created
        testthat::expect_equal("IMPUTATION_BETA" %in% colnames(imp_beta_dt), 
                               TRUE)
        #check how close it is to actual
        setkey(imp_beta_dt,SNP)
        setkey(sumstats_dt_org,SNP)
        sumstats_dt_org[imp_beta_dt,imp_beta:=i.BETA]
        #check diff as percentage of value
        sumstats_dt_org[,diff_beta:=abs((BETA-imp_beta)/BETA)]
        #expect a mean difference below 0.1%
        testthat::expect_equal(mean(sumstats_dt_org$diff_beta,na.rm = TRUE)<0.001, 
                               TRUE)
    } else {
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
        expect_equal(is_32bit_windows, TRUE)
    }
})