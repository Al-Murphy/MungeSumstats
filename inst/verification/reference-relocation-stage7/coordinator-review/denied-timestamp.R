args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[[1]])
output <- normalizePath(args[[2]])
.libPaths(c('/private/tmp/harvest-Rlib',.libPaths()))
Rcpp::sourceCpp(file.path(root,'src/preservation.cpp'))
source(file.path(root,'R/preservation.R'))
source(file.path(root,'R/projection.R'))
source(file.path(root,'tests/testthat/helper-reference-format3.R'))
library(testthat)
# Deterministically exercise the Windows-style denied-mutation branch while
# leaving all reference and normalization implementations unchanged.
Sys.setFileTime <- function(paths,time) FALSE
results <- test_file(file.path(root,'tests/testthat/test-reference-relocation.R'),
                     reporter='summary',env=.GlobalEnv)
frame <- as.data.frame(results)
write.table(frame[,c('test','nb','failed','skipped','error','warning','passed')],
            file.path(output,'denied-timestamp-summary.tsv'),sep='\t',row.names=FALSE)
cat('DENIED_TIMESTAMP blocks=',nrow(frame),' passed=',sum(frame$passed),
    ' failed=',sum(frame$failed),' errors=',sum(frame$error),
    ' warnings=',sum(frame$warning),'\n',sep='')
stopifnot(sum(frame$failed)==0,sum(frame$error)==0,sum(frame$warning)==0)
