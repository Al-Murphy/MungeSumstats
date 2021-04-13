test_that("Multi-trait GWAS handled correctly", {
  file <- tempfile()
  #write an example of a multi-trait GWAS in
  multi_trait <- c(
    "chromosome\trs_id\tmarkername\tposition_hg18\tEffect_allele\tOther_allele\tEAF_HapMapCEU\tN_SMK\tEffect_SMK\tStdErr_SMK\tP_value_SMK\tN_NONSMK\tEffect_NonSMK\tStdErr_NonSMK\tP_value_NonSMK",
    "1\t\tchr1:240307304\t240307304\tT\tG\t\t31366.8\t-0.0528\t0.0142\t0.0002004\t121792\t-0.0026\t0.0075\t0.7274",
    "1\trs1000050\tchr1:161003087\t161003087\tT\tC\t0.9\t36256.6\t0.0001\t0.0109\t0.9931\t127514\t0.0058\t0.0059\t0.3307",
    "1\trs1000073\tchr1:155522020\t155522020\tA\tG\t0.3136\t36335\t0.0046\t0.0083\t0.5812\t126780\t0.0038\t0.0045\t0.3979",
    "1\trs1000075\tchr1:94939420\t94939420\tT\tC\t0.3583\t38959.3\t-0.0013\t0.0082\t0.8687\t147567\t-0.0043\t0.0044\t0.3259",
    "1\trs1000085\tchr1:66630503\t66630503\tC\tG\t0.1667\t38761\t0.0053\t0.0095\t0.5746\t147259\t-0.0034\t0.0052\t0.5157")
  writeLines(multi_trait,con = file)
  #Run MungeSumstats code
  reformatted <- MungeSumstats::format_sumstats(file,ref_genome="GRCh37",
                                                  analysis_trait="smk")
  reformatted_res <- readLines(reformatted)
  #check manually
  multi_trait_res_smk <- c(
    "SNP\tCHR\tBP\tPOSITION_HG18\tA1\tA2\tEAF_HAPMAPCEU\tN\tBETA\tSTDERR\tP\tN_NO\tEFFECT_NO\tSTDERR_NO\tP_VALUE_NO",
    "rs1000050\t1\t161003087\t161003087\tT\tC\t0.9\t36257\t1e-04\t0.0109\t0.9931\t127514\t0.0058\t0.0059\t0.3307",
    "rs1000073\t1\t155522020\t155522020\tA\tG\t0.3136\t36335\t0.0046\t0.0083\t0.5812\t126780\t0.0038\t0.0045\t0.3979",
    "rs1000075\t1\t94939420\t94939420\tT\tC\t0.3583\t38959\t-0.0013\t0.0082\t0.8687\t147567\t-0.0043\t0.0044\t0.3259",
    "rs1000085\t1\t66630503\t66630503\tC\tG\t0.1667\t38761\t0.0053\t0.0095\t0.5746\t147259\t-0.0034\t0.0052\t0.5157"
  )
  expect_equal(multi_trait_res_smk,reformatted_res)
})
