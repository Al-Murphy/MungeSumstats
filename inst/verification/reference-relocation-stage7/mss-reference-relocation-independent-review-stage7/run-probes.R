args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(args[[1]])
evidence<-normalizePath(args[[2]])
.libPaths(c('/private/tmp/harvest-Rlib',.libPaths()))
Rcpp::sourceCpp(file.path(root,'src/preservation.cpp'))
source(file.path(root,'R/preservation.R'))
source(file.path(root,'R/projection.R'))
old<-new.env(parent=.GlobalEnv)
sys.source(file.path(evidence,'baseline-preservation.R'),old)
sys.source(file.path(evidence,'baseline-projection.R'),old)
library(testthat)
results<-test_file(file.path(evidence,'probes.R'),reporter='summary',env=.GlobalEnv)
frame<-as.data.frame(results)
write.table(frame[,c('test','nb','failed','skipped','error','warning','passed')],
  file.path(evidence,'probe-summary.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
cat('INDEPENDENT_TOTALS blocks=',nrow(frame),' assertions=',sum(frame$nb),
  ' passed=',sum(frame$passed),' failed=',sum(frame$failed),
  ' errors=',sum(frame$error),' warnings=',sum(frame$warning),'\n',sep='')
stopifnot(sum(frame$failed)==0,sum(frame$error)==0,sum(frame$warning)==0)
