cache_fixture <- function() {
    calls <- list()
    reference <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"), pos = c(10L, 20L, 30L)
    )
    data.table::setkeyv(reference, "SNP")
    loader <- function(snps, ref_genome, dbSNP, msg = NULL, chr_filt = NULL) {
        calls[[length(calls) + 1L]] <<- list(snps, ref_genome, dbSNP, chr_filt)
        data.table::copy(reference[reference$SNP %in% snps])
    }
    list(load = loader, cached = make_reference_cache(loader),
         calls = function() calls)
}

test_that("overlapping lookups cache missing IDs and isolate caller mutations", {
    f <- cache_fixture()
    x <- f$cached(c("rs2", "rs1", "rs99", "rs1"), "GRCh37", 155)
    expect_identical(x$SNP, c("rs1", "rs2"))
    data.table::set(x, j = "pos", value = 999L)
    y <- f$cached(c("rs99", "rs2", "rs1"), "GRCh37", 155)
    expect_length(f$calls(), 1L)
    expect_identical(y$pos, c(10L, 20L))
    z <- f$cached(c("rs2", "rs3", "rs99"), "GRCh37", 155)
    expect_identical(f$calls()[[2]][[1]], "rs3")
    expect_identical(z, f$load(c("rs2", "rs3", "rs99"), "GRCh37", 155))
    before <- length(f$calls())
    expect_equal(nrow(f$cached("rs99", "GRCh37", 155)), 0)
    expect_length(f$calls(), before)
})

test_that("reference builds and dbSNP versions have separate caches", {
    f <- cache_fixture()
    f$cached("rs1", "GRCh37", 155)
    f$cached("rs1", "grch37", 155)
    expect_length(f$calls(), 1L)
    f$cached("rs1", "GRCh38", 155)
    f$cached("rs1", "GRCh37", 144)
    expect_length(f$calls(), 3L)
})

test_that("empty results remain empty and are cached", {
    f <- cache_fixture()
    expect_equal(nrow(f$cached(character(), "GRCh37", 155)), 0)
    expect_equal(nrow(f$cached("rs99", "GRCh37", 155)), 0)
    before <- length(f$calls())
    expect_equal(nrow(f$cached("rs99", "GRCh37", 155)), 0)
    expect_length(f$calls(), before)
})

test_that("unusual IDs and chromosome restrictions use the original loader", {
    f <- cache_fixture()
    f$cached("rs01", "GRCh37", 155)
    f$cached(NA_character_, "GRCh37", 155)
    f$cached("rs1", "GRCh37", 155, chr_filt = "1")
    f$cached("rs1", "GRCh37", c(144, 155))
    f$cached("rs1", "GRCh37", NA_real_)
    expect_length(f$calls(), 5L)
    expect_identical(f$calls()[[3]][[4]], "1")
})

test_that("alias expansion falls back without caching an unsafe result", {
    calls <- 0L
    loader <- function(snps, ...) {
        calls <<- calls + 1L
        data.table::data.table(SNP = "rs2", pos = 20L)
    }
    cached <- make_reference_cache(loader)
    expect_identical(cached("rs1", "GRCh37", 155), loader("rs1"))
    before <- calls
    cached("rs1", "GRCh37", 155)
    expect_equal(calls - before, 2L)
})

test_that("format_sumstats restores cache state after nested calls and errors", {
    # Isolate the function's bindings without loading reference packages.
    scope <- new.env(parent = environment(format_sumstats))
    scope$.reference_cache <- new.env(parent = emptyenv())
    scope$.reference_cache$loader <- NULL
    scope$format_sumstats <- format_sumstats
    environment(scope$format_sumstats) <- scope
    outer <- NULL
    nested <- FALSE
    scope$check_save_path <- function(...) {
        if (!nested) {
            outer <<- scope$.reference_cache$loader
            expect_true(is.function(outer))
            nested <<- TRUE
            expect_error(scope$format_sumstats("unused", cache_reference = FALSE),
                         "nested exit")
            expect_identical(scope$.reference_cache$loader, outer)
            stop("outer exit")
        }
        expect_null(scope$.reference_cache$loader)
        stop("nested exit")
    }
    expect_error(scope$format_sumstats("unused", cache_reference = TRUE),
                 "outer exit")
    expect_null(scope$.reference_cache$loader)
    expect_error(scope$format_sumstats("unused", cache_reference = NA),
                 "cache_reference must be TRUE or FALSE")
    expect_false(formals(format_sumstats)$cache_reference)
})

test_that("the lookup entry point dispatches to the active cache only", {
    scope <- new.env(parent = environment(load_ref_genome_data))
    scope$.reference_cache <- new.env(parent = emptyenv())
    scope$.reference_cache$loader <- NULL
    scope$load_ref_genome_data_uncached <- function(...) "uncached"
    lookup <- load_ref_genome_data
    environment(lookup) <- scope
    expect_identical(lookup("rs1", "GRCh37", 155), "uncached")
    scope$.reference_cache$loader <- function(...) "cached"
    expect_identical(lookup("rs1", "GRCh37", 155), "cached")
    scope$.reference_cache$loader <- NULL
    expect_identical(lookup("rs1", "GRCh37", 155), "uncached")
})
