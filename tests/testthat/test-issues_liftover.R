# Regression tests for GitHub issues on genome-build liftover (offline: the
# bundled Ensembl chain is used). Issues: #84, #104, #107, #118, #128, #133,
# #136, #183, #184, #197, #211.
# A fresh copy of a bundled chain in its own temp dir, so liftover() never
# touches the installed file and no tempdir() cache is involved.
fresh_chain <- function(f = "GRCh37_to_GRCh38.chain.gz", gunzip = FALSE) {
    d <- tempfile("chain_")
    dir.create(d)
    p <- file.path(d, f)
    file.copy(system.file("extdata", f, package = "MungeSumstats"), p)
    if (gunzip) p <- as.character(R.utils::gunzip(p))
    stopifnot(file.exists(p))
    p
}
# eduAttainOkbay with standard headers (cheaper than formatted_example())
edu_dt <- function() {
    dt <- data.table::fread(
        system.file("extdata", "eduAttainOkbay.txt", package = "MungeSumstats"),
        nThread = 1
    )
    data.table::setnames(dt, c("MarkerName", "POS", "EAF", "Beta", "Pval"),
                         c("SNP", "BP", "FRQ", "BETA", "P"))
    dt[, CHR := as.character(CHR)]
    dt
}
# Turn any chain download into an error: get_chain_file() announces each
# download through messager() before it calls download.file()
no_downloads <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(
        messager = function(..., v = TRUE) {
            msg <- paste(...)
            if (grepl("Downloading", msg)) stop("Download attempted: ", msg)
            if (v) message(msg)
        },
        .package = "MungeSumstats", .env = env
    )
}
# run liftover() and keep both its result and its messages
lift_msgs <- function(...) {
    out <- NULL
    msgs <- testthat::capture_messages(out <- MungeSumstats::liftover(...))
    list(out = out, msgs = msgs)
}

test_that("#104: a chrX SNP lifts to the right CHR and position", {
    no_downloads()
    dt <- data.table::data.table(
        SNP = c("rs190850470", "rs12987662"), CHR = c("X", "2"),
        BP = c(50858748L, 100821548L), A1 = c("G", "A"), A2 = c("A", "C"),
        BETA = c(-0.0303521, 0.027), SE = c(0.0670233, 0.003),
        P = c(0.650651, 2.693e-24)
    )
    res <- lift_msgs(sumstats_dt = data.table::copy(dt), ref_genome = "GRCh37",
                     convert_ref_genome = "GRCh38", imputation_ind = FALSE,
                     local_chain = fresh_chain())
    expect_false(any(grepl("Downloading", res$msgs)))
    out <- res$out
    expect_equal(nrow(out), 2)
    expect_equal(as.character(out[SNP == "rs190850470"]$CHR), "X")
    expect_equal(out[SNP == "rs190850470"]$BP, 51115913L)
    expect_equal(as.character(out[SNP == "rs12987662"]$CHR), "2")
    expect_equal(out[SNP == "rs12987662"]$BP, 100205086L)
})

test_that("#107/#118/#136: no conversion requested -> unchanged and silent", {
    dt <- data.table::data.table(SNP = c("rs1", "rs2"), CHR = c("1", "2"),
                                 BP = c(100L, 200L), A1 = "A", A2 = "G",
                                 P = 0.1)
    expect_silent(out <- MungeSumstats::liftover(
        sumstats_dt = data.table::copy(dt), ref_genome = "GRCh37",
        convert_ref_genome = NULL
    ))
    expect_equal(out, dt)
    # same build: no "Could not recognize genome build" / liftover
    expect_silent(out <- MungeSumstats::liftover(
        sumstats_dt = data.table::copy(dt), ref_genome = "GRCh37",
        convert_ref_genome = "hg19"
    ))
    expect_equal(out, dt)
    # #118/#136: liftover is exported
    expect_true("liftover" %in% getNamespaceExports("MungeSumstats"))
})

test_that("#128/#133: Ensembl is the default chain source; no UCSC shipped", {
    expect_equal(formals(MungeSumstats::format_sumstats)$chain_source,
                 "ensembl")
    expect_equal(formals(MungeSumstats::liftover)$chain_source, "ensembl")
    shipped <- list.files(system.file("extdata", package = "MungeSumstats"),
                          pattern = "chain")
    expect_setequal(shipped, c("GRCh37_to_GRCh38.chain.gz",
                               "GRCh38_to_GRCh37.chain.gz"))
    expect_equal(
        MungeSumstats:::.get_ensembl_chain_remote("hg19", "hg38"),
        paste0("ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/",
               "GRCh37_to_GRCh38.chain.gz")
    )
    expect_equal(
        MungeSumstats:::.get_ucsc_chain_remote("hg38", "hg19"),
        paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/",
               "hg38ToHg19.over.chain.gz")
    )
    expect_error(MungeSumstats:::get_chain_file(from = "hg19", to = "hg19"),
                 "Cannot get a chain file from one reference to the same")
})

test_that("#84/#128: get_chain_file uses an existing chain, no download", {
    no_downloads()
    ens <- fresh_chain(gunzip = TRUE)
    # see the known issue below: a temp path with "<char>gz" in it would
    # hide the cached file
    skip_if(grepl(".gz", dirname(ens)), "temp path matches '.gz'")
    msgs <- testthat::capture_messages(
        ch <- MungeSumstats:::get_chain_file(from = "hg19", to = "hg38",
                                             chain_source = "ensembl",
                                             save_dir = dirname(ens))
    )
    expect_true(any(grepl("Using existing chain file", msgs)))
    expect_false(any(grepl("Downloading", msgs)))
    expect_s4_class(ch, "Chain")
    expect_gt(length(ch), 0)
    # UCSC chains carry a commercial-use licence warning
    ucsc_dir <- tempfile("ucsc_")
    dir.create(ucsc_dir)
    stopifnot(file.copy(ens, file.path(ucsc_dir, "hg19ToHg38.over.chain")))
    msgs <- testthat::capture_messages(
        MungeSumstats:::get_chain_file(from = "hg19", to = "hg38",
                                       chain_source = "ucsc",
                                       save_dir = ucsc_dir)
    )
    expect_true(any(grepl("licence for commercial use", msgs)))
    expect_false(any(grepl("Downloading", msgs)))
})

test_that("#84: get_chain_file finds a cached chain in any save_dir", {
    no_downloads()
    save_dir <- file.path(tempfile("chain_"), "gwas_gz_files")
    dir.create(save_dir, recursive = TRUE)
    stopifnot(file.copy(fresh_chain(gunzip = TRUE),
                        file.path(save_dir, "GRCh37_to_GRCh38.chain")))
    expect_message(
        ch <- MungeSumstats:::get_chain_file(from = "hg19", to = "hg38",
                                             chain_source = "ensembl",
                                             save_dir = save_dir),
        "Using existing chain file"
    )
    expect_s4_class(ch, "Chain")
})

test_that("#183/#184/#197: local_chain .gz is kept, reusable; extra cols", {
    no_downloads()
    dt <- edu_dt()
    dt[, ROW_NUMBER := .I]
    dt[, ANNOTATION := paste0("a", .I)]
    chain_gz <- fresh_chain()
    lift <- function() {
        lift_msgs(sumstats_dt = data.table::copy(dt), ref_genome = "hg19",
                  convert_ref_genome = "hg38", local_chain = chain_gz)
    }
    res1 <- lift()
    expect_true(any(grepl("Using local chain file", res1$msgs)))
    expect_false(any(grepl("Downloading", res1$msgs)))
    out1 <- res1$out
    expect_equal(nrow(out1), nrow(dt))
    expect_setequal(out1$SNP, dt$SNP)
    expect_true(all(out1[order(SNP)]$BP != dt[order(SNP)]$BP))
    expect_true("IMPUTATION_gen_build" %in% names(out1))
    # #197: user columns travel with their rows
    expect_equal(out1[order(SNP)]$ROW_NUMBER, dt[order(SNP)]$ROW_NUMBER)
    expect_equal(out1[order(SNP)]$ANNOTATION, dt[order(SNP)]$ANNOTATION)
    # #183: the user's .gz survives and can be used again
    expect_true(file.exists(chain_gz))
    expect_identical(lift()$out, out1)
    expect_error(MungeSumstats::liftover(sumstats_dt = data.table::copy(dt),
                                         ref_genome = "hg19",
                                         convert_ref_genome = "hg38",
                                         local_chain = tempfile()),
                 "local_chain parameter is invalid")
})

test_that("#211: CHR 23 is lifted as X and rmv_chr removes rows", {
    no_downloads()
    dt <- edu_dt()
    # rs301800 (1:8490603) relabelled as an X SNP inside a chain block
    dt[SNP == "rs301800",
       `:=`(SNP = "rs796052235", CHR = "23", A1 = "G", A2 = "A")]
    out <- suppressMessages(MungeSumstats::liftover(
        sumstats_dt = data.table::copy(dt), ref_genome = "hg19",
        convert_ref_genome = "hg38", local_chain = fresh_chain()
    ))
    expect_equal(nrow(out), nrow(dt))
    expect_false("23" %in% as.character(out$CHR))
    expect_equal(as.character(out[SNP == "rs796052235"]$CHR), "X")
    # X:7673882 -> X:7705841 block: 8490603 - 7673882 + 7705841
    expect_equal(out[SNP == "rs796052235"]$BP, 8522562L)
    expect_message(
        out <- MungeSumstats::liftover(
            sumstats_dt = data.table::copy(dt), ref_genome = "hg19",
            convert_ref_genome = "hg38", local_chain = fresh_chain(),
            rmv_chr = "X"
        ),
        "1 SNPs are on chromosomes X and will be removed"
    )
    expect_equal(nrow(out), nrow(dt) - 1)
    expect_false("X" %in% as.character(out$CHR))
})
