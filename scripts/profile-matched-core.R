#!/usr/bin/env Rscript
# Same arguments as benchmark-matched.R, plus an eighth runner snapshot path.
args <- commandArgs(TRUE)
if(length(args)!=8) stop('mode input reference output metrics cache threads runner_snapshot')
runner <- args[[8]]
stage <- '/private/tmp/gwas-harvest-stage2/Rlib'
.libPaths(c(if(args[[1]]=='legacy') file.path(stage,'legacy') else Sys.getenv('MSS_NEW_RLIB',stage),stage,'/private/tmp/harvest-Rlib',.libPaths()))
library(MungeSumstats)
.mss_profile <- new.env(parent=emptyenv())
.mss_profile$times <- list()
targets <- if(args[[1]]=='legacy') c('check_allele_flip','check_bi_allelic') else 'mss_harmonise_cpp'
ns <- asNamespace('MungeSumstats')
for(name in targets) {
    label <- name
    # Wrapping the real boundary changes no implementation or scientific choice.
    wrapper <- local({original <- get(name,ns);key <- label
        function(...) {
            start <- proc.time()
            on.exit({delta<-proc.time()-start
                     previous<-.mss_profile$times[[key]]
                     if(is.null(previous)) previous<-c(user.self=0,sys.self=0,elapsed=0)
                     .mss_profile$times[[key]]<-previous+delta[names(previous)]},add=TRUE)
            original(...)
        }
    })
    unlockBinding(name,ns);assign(name,wrapper,ns);lockBinding(name,ns)
}
# The runner accepts only seven args; evaluate its immutable text after replacing
# commandArgs in an isolated environment, preserving the exact actual argv.
environment <- new.env(parent=globalenv())
environment$commandArgs <- function(trailingOnly=FALSE) if(trailingOnly) args[1:7] else base::commandArgs(FALSE)
sys.source(runner,envir=environment)
core <- do.call(rbind,lapply(names(.mss_profile$times),function(name)
    data.frame(boundary=name,t(.mss_profile$times[[name]]),row.names=NULL)))
write.table(core,paste0(args[[5]],'.core.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
cat('Instrumented core boundaries; warm totals include both priming and processing.\n')
