# Small end-to-end benchmark-CLI validation; no performance inference.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[[1]]) else normalizePath('.')
test_lib <- Sys.getenv('MSS_TEST_RLIB')
if(nzchar(test_lib)) .libPaths(c(test_lib,.libPaths()))
Rcpp::sourceCpp(file.path(root,'src/preservation.cpp'))
source(file.path(root,'R/preservation.R'))
d <- tempfile('mss-benchmark-test-');dir.create(d)
writeLines(c('>1','GCGT'),file.path(d,'reference.fa'))
writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t1\tG\tA,T\trs1\told1'),file.path(d,'reference.tsv'))
reference <- prepare_reference(file.path(d,'reference.fa'),file.path(d,'reference.tsv'),file.path(d,'reference'),build='GRCh38',version='benchmark-fixture',batch_size=1)
writeLines(c('CHR\tBP\tSNP\tA1\tA2\tBETA\tFRQ','1\t1\told1\tG\tT\t.1\t0','1\t1\trs1\tA\tT\t.3\t.2'),file.path(d,'input.tsv'))
run_cli <- function(rows,edition,metrics) {
  argv <- c(file.path(root,'scripts/benchmark-preservation.R'),
            '--reference',reference$path,'--input',file.path(d,'input.tsv'),
            '--output',file.path(d,edition),'--metrics',file.path(d,metrics),
            '--rows',rows,'--batch-size','1','--frequency-type','EAF',
            '--effect-type','beta','--cache-state','uncontrolled-smoke',
            '--source-root',root)
  suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),shQuote(argv),stdout=TRUE,stderr=TRUE))
}
output <- run_cli('2','edition','metrics.tsv')
if(!is.null(attr(output,'status')))stop(paste(output,collapse='\n'))
metrics <- read.delim(file.path(d,'metrics.tsv'),colClasses='character')
value <- function(key) metrics$value[metrics$key==key]
stopifnot(identical(value('source_rows'),'2'),identical(value('expected_rows'),'2'),
          identical(value('count_validation'),'pass'),identical(value('rss_measured'),'false'),
          identical(value('cache_state'),'uncontrolled-smoke'),
          nzchar(value('input_checksum')),as.numeric(value('harmonise_elapsed_seconds'))>=0)
failed <- run_cli('3','count-mismatch','mismatch.tsv')
stopifnot(!is.null(attr(failed,'status')))
mismatch <- read.delim(file.path(d,'mismatch.tsv'),colClasses='character')
stopifnot(identical(mismatch$value[mismatch$key=='count_validation'],'fail'))
unlink(d,recursive=TRUE)
cat('benchmark CLI smoke fixtures passed\n')
