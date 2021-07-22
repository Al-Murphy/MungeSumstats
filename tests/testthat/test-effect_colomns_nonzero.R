test_that("Filter effect columns so !=0", {
    file <- tempfile()
    #write the Educational Attainment GWAS to a temp file for testing
    eduAttainOkbay <- readLines(system.file("extdata","eduAttainOkbay.txt",
                                            package="MungeSumstats"))
    writeLines(eduAttainOkbay,con = file)
    #read it in and make N
    sumstats_dt <- data.table::fread(file, nThread = 1)
    #Change BETA to 0 for some SNPs
    num_snps <- 5
    sumstats_dt[seq_len(num_snps),Beta:=0]
    #get SNPs names
    rmv_snps <- sumstats_dt[Beta==0,]$MarkerName
    data.table::fwrite(x=sumstats_dt, file=file, sep="\t")
    #Run MungeSumstats code
    reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                  effect_columns_nonzero=TRUE,
                                                  on_ref_genome = FALSE,
                                                  strand_ambig_filter=FALSE,
                                                  bi_allelic_filter=FALSE,
                                                  allele_flip_check=FALSE)
    res_dt <- data.table::fread(reformatted, nThread = 1)
    testthat::expect_equal(all(!rmv_snps %in% res_dt$SNP),TRUE)
})
