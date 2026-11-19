# Verify generated .Call registration without loading the full package.
test_lib <- Sys.getenv("MSS_TEST_RLIB")
if (nzchar(test_lib)) .libPaths(c(test_lib, .libPaths()))
library(Rcpp)
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[[1]]) else normalizePath(".")
build <- tempfile('mss-shlib-');dir.create(build)
files <- c('preservation.cpp','preservation_core.h','preservation_stream.h','projection.h','RcppExports.cpp')
file.copy(file.path(root,'src',files),build)
Sys.setenv(PKG_CPPFLAGS=paste0('-I',shQuote(system.file('include',package='Rcpp'))),PKG_CXXFLAGS='-std=c++17')
old <- setwd(build)
status <- system2(file.path(R.home('bin'),'R'),c('CMD','SHLIB','preservation.cpp','RcppExports.cpp','-o',paste0('MungeSumstats',.Platform$dynlib.ext)))
stopifnot(status==0)
dll <- dyn.load(file.path(build,paste0('MungeSumstats',.Platform$dynlib.ext)))
expected <- c(
    '_MungeSumstats_mss_prepare_cpp' = 5L,
    '_MungeSumstats_mss_open_reference_cpp' = 2L,
    '_MungeSumstats_mss_harmonise_cpp' = 15L,
    '_MungeSumstats_mss_parser_open_cpp' = 2L,
    '_MungeSumstats_mss_parser_feed_cpp' = 3L,
    '_MungeSumstats_mss_parser_batch_cpp' = 8L,
    '_MungeSumstats_mss_write_normalized_cpp' = 5L,
    '_MungeSumstats_mss_projection_open_cpp' = 5L,
    '_MungeSumstats_mss_project_cpp' = 10L
)
registered <- getDLLRegisteredRoutines(dll)$.Call
stopifnot(
    'Unexpected registered .Call routine names' =
        identical(sort(names(registered)), sort(names(expected))),
    'Unexpected registered .Call routine arities' =
        identical(vapply(registered[names(expected)], function(routine) routine$numParameters,
                         integer(1)), expected)
)
for (name in names(registered)) {
    assign(name, getNativeSymbolInfo(name, dll))
}
source(file.path(root,'R','RcppExports.R'))
source(file.path(root,'R','preservation.R'))
source(file.path(root,'R','projection.R'))
library(testthat)
source(file.path(root,'tests/testthat/test-preservation.R'))
source(file.path(root,'tests/testthat/test-source-reference-contract.R'))
source(file.path(root,'tests/testthat/test-equivalent-indel-context.R'))
source(file.path(root,'tests/testthat/test-quoted-stream.R'))
source(file.path(root,'tests/testthat/test-preservation-stream-performance.R'))
source(file.path(root,'tests/testthat/test-projection.R'))
source(file.path(root,'tests/testthat/helper-reference-format3.R'))
source(file.path(root,'tests/testthat/test-reference-format3.R'))
source(file.path(root,'tests/testthat/test-reference-relocation.R'))
source(file.path(root,'tests/testthat/test-fai-safety.R'))
setwd(old)
cat('registered generated bindings fixtures passed\n')
