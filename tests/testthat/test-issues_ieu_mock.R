# Regression tests for GitHub issues on OpenGWAS import and dbSNP tarballs.
# Fully offline (ieugwasr/downloader/BSgenome calls are mocked).
# Issues: #52, #66, #74, #109, #204, #209, #210.
test_that("#52/#109: find_sumstats computes N per row", {
    meta <- data.frame(
        id = c("a", "b", "c", "d"), trait = c("T1", "T2", "T3", "T4"),
        year = 2020, consortium = "C", author = "A", population = "European",
        category = "Disease", subcategory = "Immune", build = "HG19/GRCh37",
        pmid = 1, nsnp = 1e6, sample_size = c(NA, 500, NA, NA),
        ncase = c(100, NA, NA, 1000), ncontrol = c(200, NA, NA, 2000)
    )
    testthat::local_mocked_bindings(
        gwasinfo = function(id = NULL, ...) {
            if (is.null(id)) meta else meta[meta$id %in% id, ]
        },
        .package = "ieugwasr"
    )
    res <- suppressMessages(MungeSumstats::find_sumstats())
    # N = ncase + ncontrol per row (#109: not summed over all rows)
    expect_equal(res[id == "a"]$N, 300)
    expect_equal(res[id == "d"]$N, 3000)
    # sample_size takes precedence; no counts at all -> NA, not 0
    expect_equal(res[id == "b"]$N, 500)
    expect_true(is.na(res[id == "c"]$N))
    res <- suppressMessages(MungeSumstats::find_sumstats(ids = c("a", "b")))
    expect_setequal(res$id, c("a", "b"))
    expect_true("query" %in% names(res))
})

test_that("#66/#74: import_sumstats reuses an existing formatted file", {
    id <- "ieu-a-298"
    testthat::local_mocked_bindings(
        gwasinfo_files = function(id, ...) {
            list(paste0("https://example.org/files/", id, "/", id,
                        c(".vcf.gz", ".vcf.gz.tbi")))
        },
        .package = "ieugwasr"
    )
    testthat::local_mocked_bindings(
        download_vcf = function(...) stop("should not download"),
        .package = "MungeSumstats"
    )
    save_dir <- tempfile("ieu_")
    sp <- file.path(save_dir, id, paste0(id, ".tsv.gz"))
    dir.create(dirname(sp), recursive = TRUE)
    dt <- data.table::data.table(SNP = c("rs1", "rs2"), CHR = 1L,
                                 BP = c(100L, 200L), A1 = "A", A2 = "G",
                                 BETA = 0.1, SE = 0.01, P = 0.05)
    suppressMessages(MungeSumstats::write_sumstats(dt, save_path = sp))
    # #74: no download and no re-munging when <save_dir>/<id>/<id>.tsv.gz
    # exists. This path relies on format_sumstats() never forcing `path`
    # (vcf_paths is never assigned); a failure there would come back as an
    # error object, not be thrown.
    expect_message(
        res <- MungeSumstats::import_sumstats(ids = id, save_dir = save_dir),
        "Importing previously formatted file"
    )
    expect_false(methods::is(res[[id]], "error"))
    expect_equal(res[[id]], sp)
    # #66: one directory per dataset, with its own logs folder
    expect_true(dir.exists(file.path(save_dir, id, "logs")))
})

test_that("#209/#210: download_vcf fetches the index from its own URL", {
    calls <- character()
    testthat::local_mocked_bindings(
        downloader = function(input_url, output_path, ...) {
            calls <<- c(calls, input_url)
            f <- file.path(output_path, basename(input_url))
            file.create(f)
            f
        },
        .package = "MungeSumstats"
    )
    vcf_dir <- tempfile("vcf_")
    dir.create(vcf_dir)
    # OpenGWAS signs each file URL separately, so the .tbi URL is not
    # simply paste0(vcf_url, ".tbi")
    urls <- c("https://example.org/tok1/ieu-a-298.vcf.gz",
              "https://example.org/tok2/ieu-a-298.vcf.gz.tbi")
    out <- suppressMessages(MungeSumstats::download_vcf(
        vcf_urls = urls, vcf_dir = vcf_dir, force_new = TRUE
    ))
    expect_equal(calls, urls)
    expect_equal(out$save_path, file.path(vcf_dir, "ieu-a-298.vcf.gz"))
    expect_equal(out$index_path, file.path(vcf_dir, "ieu-a-298.vcf.gz.tbi"))
    # second call reuses the download
    expect_message(
        out2 <- MungeSumstats::download_vcf(vcf_urls = urls,
                                            vcf_dir = vcf_dir),
        "Using previously downloaded VCF"
    )
    expect_equal(out2, out)
    expect_length(calls, 2)
    # a local VCF is used in place
    local_vcf <- tempfile(fileext = ".vcf.gz")
    file.copy(system.file("extdata", "ALSvcf.vcf.gz",
                          package = "MungeSumstats"), local_vcf)
    expect_message(out3 <- MungeSumstats::download_vcf(vcf_urls = local_vcf),
                   "Using local VCF")
    expect_equal(out3$save_path, local_vcf)
    expect_length(calls, 2)
})

test_that("#204: each dbSNP tarball gets its own extraction", {
    # BSgenome::OnDiskLongTable() would need real SNPlocs data: return the
    # extdata directory it was given instead
    testthat::local_mocked_bindings(OnDiskLongTable = function(...) ..1,
                                    .package = "BSgenome")
    src <- tempfile("tarballs_")
    dir.create(src)
    make_tarball <- function(build) {
        pkg <- paste0("SNPlocs.Hsapiens.dbSNP156.", build)
        dir.create(file.path(src, pkg, "inst", "extdata"), recursive = TRUE)
        writeLines(build, file.path(src, pkg, "inst", "extdata", "id.txt"))
        # unique name so no extraction from another test can be reused
        tb <- file.path(src, paste0(basename(tempfile(paste0(pkg, "_"))),
                                    ".tar.gz"))
        owd <- setwd(src)
        on.exit(setwd(owd))
        utils::tar(tb, files = pkg, compression = "gzip")
        tb
    }
    tb37 <- make_tarball("GRCh37")
    tb38 <- make_tarball("GRCh38")
    load_tb <- function(build, tb) {
        MungeSumstats::load_snp_loc_data(ref_genome = build, dbSNP = 156,
                                         dbSNP_tarball = tb)
    }
    expect_message(p37 <- load_tb("GRCh37", tb37),
                   "Loading SNPlocs data from tarball")
    p38 <- suppressMessages(load_tb("GRCh38", tb38))
    expect_equal(readLines(file.path(p37, "id.txt")), "GRCh37")
    # a second tarball in the same session must not reuse the first one
    expect_equal(readLines(file.path(p38, "id.txt")), "GRCh38")
    # the same tarball is only extracted once per session
    expect_equal(suppressMessages(load_tb("GRCh37", tb37)), p37)
    unlink(c(src, file.path(tempdir(), paste0("snp_loc_", basename(
        c(tb37, tb38))))), recursive = TRUE)
})
