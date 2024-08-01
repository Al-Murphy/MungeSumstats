test_that("Test infer effect column function works", {
  ## Call uses reference genome as default with more than 2GB of memory,
  ## which is more than what 32-bit Windows can handle so remove tests
  is_32bit_windows <-
    .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
  if (!is_32bit_windows) {
    a <- MungeSumstats::formatted_example()
    b <- copy(a)
    
    #check if inferring effect direction is working
    
    #make A1 the effect col not A2
    #change FRQ and BETA to reflect this change
    data.table::setnames(b,"A2","A3")
    data.table::setnames(b,"A1","A2")
    data.table::setnames(b,"A3","A1")
    c <- copy(b)
    data.table::setnames(b,"BETA","BETA1")
    data.table::setnames(c,"FRQ","A1FRQ")
    
    # if the allele columns are named unambiguously, formatting will be performed 
    # correctly - use this to test against as eff will be flipped compared to org
    b_renamed <- copy(b)
    data.table::setnames(b_renamed, "A1", "effect_allele")
    data.table::setnames(b_renamed, "A2", "non_effect_allele")
    b_renamed_for <- MungeSumstats::format_sumstats(b_renamed, return_data = TRUE,
                                                    #all just make MSS run faster
                                                    ref_genome = 'GRCh37',
                                                    on_ref_genome = FALSE,
                                                    strand_ambig_filter = FALSE,
                                                    bi_allelic_filter = FALSE,
                                                    allele_flip_check = FALSE,
                                                    dbSNP=144)
    
    # what if BETA gives the direction
    b_for <- MungeSumstats::format_sumstats(b, return_data = TRUE, 
                                            #all just make MSS run faster
                                            ref_genome = 'GRCh37',
                                            on_ref_genome = FALSE,
                                            strand_ambig_filter = FALSE,
                                            bi_allelic_filter = FALSE,
                                            allele_flip_check = FALSE,
                                            dbSNP=144)
    
    # what if FRQ gives the direction
    c_for <- MungeSumstats::format_sumstats(c, return_data = TRUE, 
                                            #all just make MSS run faster
                                            ref_genome = 'GRCh37',
                                            on_ref_genome = FALSE,
                                            strand_ambig_filter = FALSE,
                                            bi_allelic_filter = FALSE,
                                            allele_flip_check = FALSE,
                                            dbSNP=144)
    
    testthat::expect_equal(
      all.equal(b_renamed_for, b_for,ignore.col.order = TRUE),
      TRUE
    )
    
    testthat::expect_equal(
      all.equal(b_renamed_for, c_for,ignore.col.order = TRUE),
      TRUE
    )
    
    
    #finally check if the ref genome can be used to infer rather than being 
    #told in eff/frq cols - just subset for speed
    snps <- c("rs11210860","rs34305371","rs1008078","rs11588857","rs1777827",
              "rs76076331","rs2457660","rs10496091","rs4851251","rs12987662",
              "rs10930008","rs301800","rs2568955","rs61787263","rs2992632",
              "rs11689269","rs11690172")
    d <- copy(b)
    data.table::setnames(d,"BETA1","BETA")
    d_for <- MungeSumstats::format_sumstats(d[!SNP %in% snps,], return_data = TRUE, 
                                            on_ref_genome = TRUE,
                                            #all just make MSS run faster
                                            ref_genome = 'GRCh37',
                                            strand_ambig_filter = FALSE,
                                            bi_allelic_filter = FALSE,
                                            allele_flip_check = FALSE,
                                            dbSNP=144)
    data.table::setkey(b_renamed_for,"SNP")
    data.table::setkey(d_for,"SNP")
    testthat::expect_equal(
      all.equal(b_renamed_for[!SNP %in% snps], d_for,
                ignore.col.order = TRUE),
      TRUE
    )
  }    
  else{
    testthat::expect_equal(is_32bit_windows, TRUE)
    testthat::expect_equal(is_32bit_windows, TRUE)
    testthat::expect_equal(is_32bit_windows, TRUE)
  }
})
