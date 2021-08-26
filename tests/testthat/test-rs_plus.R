test_that("Handle extra info in rs id column", {
    file <- tempfile()
    #create dataset to test
    set.seed(101)
    a <- data.table::data.table(SNP=c("rs140052487:C:A","rs558796213:G:T",
                                      "rs561234294:A:G","rs2462492:C:T",
                                      "rs548455890:T:G"),
                                "P"=(runif(5)/10),
                                "FRQ"=runif(5),"BETA"=runif(5),"CHR"=rep(1,5),
                                "BP"=c(54353,54564,54591,54676,54763))
    #Run MungeSumstats code
    reformatted <- MungeSumstats::format_sumstats(a,ref_genome="GRCh37",
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter=FALSE,
                                                  bi_allelic_filter=FALSE,
                                                  allele_flip_check=FALSE,
                                                  log_folder_ind = TRUE,
                                                  imputation_ind = TRUE,
                                                  return_data = TRUE)
    
    #Should equal separated version
    set.seed(101)
    b <- data.table::data.table(SNP=c("rs140052487","rs558796213","rs561234294",
                          "rs2462492","rs548455890"),"P"=(runif(5)/10),
                    "FRQ"=runif(5),"BETA"=runif(5),"CHR"=rep(1,5),
                    "BP"=c(54353,54564,54591,54676,54763),
                    "A1"=c("C","G","A","C","T"),"A2"=c("A","T","G","T","G")
                    )
    reformatted2 <- MungeSumstats::format_sumstats(b,ref_genome="GRCh37",
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter=FALSE,
                                                  bi_allelic_filter=FALSE,
                                                  allele_flip_check=FALSE,
                                                  log_folder_ind = TRUE,
                                                  imputation_ind = TRUE,
                                                  return_data = TRUE)
    #should be the same
    expect_equal(setequal(reformatted$sumstats,reformatted2$sumstats),TRUE)
})
