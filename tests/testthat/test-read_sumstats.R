test_that("Test that different kinds of files can be read in.", {
    ### read data
    path <- system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats")
    sumstats_dt <- data.table::fread(path, nThread = 1)
    sumstats_dt <- standardise_sumstats_column_headers_crossplatform(
        sumstats_dt = sumstats_dt,
        mapping_file = sumstatsColHeaders
    )[["sumstats_dt"]]
    col_order <- colnames(sumstats_dt)
    ### Test tsv
    tsv_ss <- sumstats_dt
    data.table::setkey(tsv_ss, SNP)
    ### Test csv
    path_csv <- tempfile()
    data.table::fwrite(sumstats_dt, path_csv, sep = ",")
    csv_ss <- MungeSumstats::read_sumstats(path = path_csv)
    data.table::setkey(csv_ss, SNP)
    #### Test space-separated
    path_space <- tempfile()
    data.table::fwrite(sumstats_dt, path_space, sep = " ")
    space_ss <- MungeSumstats::read_sumstats(path = path_space)
    data.table::setkey(space_ss, SNP)
    ### Test VCF
    vcf_tmp <- tempfile(fileext = ".vcf.gz")
    write_sumstats(
        sumstats_dt = sumstats_dt,
        save_path = vcf_tmp, sep = "\t", write_vcf = TRUE
    )
    vcf_ss <- MungeSumstats::read_sumstats(path = vcf_tmp)
    vcf_ss <- standardise_sumstats_column_headers_crossplatform(
        sumstats_dt = vcf_ss,
        mapping_file = sumstatsColHeaders
    )[["sumstats_dt"]]
    empty_cols <- names(check_empty_cols(sumstats_file = vcf_ss))
    #vcf_ss[, (empty_cols) := NULL]
    vcf_ss[, INFO := NULL]
    for (x in c("FRQ", "BETA", "SE", "P")) {
        vcf_ss[, (x) := as.numeric(get(x))]
    }
    data.table::setcolorder(vcf_ss, col_order)
    data.table::setkey(vcf_ss, SNP)

    #### Test that all dts are the same
    objs <- mget(c("tsv_ss", "csv_ss", "space_ss", "vcf_ss"))
    iterative_tests <- outer(objs, objs, Vectorize(all.equal))
    if (any(is(unlist(iterative_tests[1, ]), "character"))) {
        all_same <- FALSE
    } else {
        all_same <- all(colSums(iterative_tests) == nrow(iterative_tests))
    }

    expect_equal(all_same, TRUE)
})
