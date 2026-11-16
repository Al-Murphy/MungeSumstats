script_start <- proc.time()
a <- commandArgs(TRUE)
if(length(a)!=1L)stop('One immutable run specification is required')
.libPaths(c('/private/tmp/harvest-Rlib',.libPaths()))
spec <- jsonlite::read_json(a[[1]],simplifyVector=TRUE)
.libPaths(c(spec$library,'/private/tmp/gwas-harvest-stage2/Rlib/reference-contigs-v1',
            '/private/tmp/gwas-harvest-stage2/Rlib','/private/tmp/harvest-Rlib',.libPaths()))
library(Rcpp)
suppressPackageStartupMessages(library(MungeSumstats))
library(data.table)
setDTthreads(as.integer(spec$threads))
if(normalizePath(find.package('MungeSumstats'))!=normalizePath(file.path(spec$library,'MungeSumstats')))
    stop('Unexpected installed package selected')
dll <- getLoadedDLLs()[['MungeSumstats']][['path']]
if(normalizePath(dll)!=normalizePath(spec$dll_path)||unname(tools::md5sum(dll))!=spec$dll_md5)
    stop('Loaded DLL differs from frozen profile')
if(file.exists(spec$output)||file.exists(spec$metrics))stop('Existing output refused')
metrics <- list(schema=1L,version=spec$version,expected_implementation=spec$implementation,
    source_production_revision=spec$production_revision,library=normalizePath(spec$library),
    loaded_dll_path=normalizePath(dll),loaded_dll_md5=unname(tools::md5sum(dll)),
    loaded_dll_sha256_from_verified_profile=spec$dll_sha256,
    input=normalizePath(spec$input),reference=normalizePath(spec$reference),
    package_version=as.character(packageVersion('MungeSumstats')),R_version=R.version.string,
    threads=spec$threads,batch_rows=spec$batch_rows,projection='not requested',
    workload=spec$workload,cache_caveat=spec$cache_caveat)
setup <- proc.time()-script_start
for(k in c('elapsed','user.self','sys.self'))metrics[[paste0('package_setup_',k)]]<-unname(setup[[k]])
stamp <- function(name,expr) {
    time <- system.time(value <- force(expr))
    for(k in c('elapsed','user.self','sys.self'))metrics[[paste(name,k,sep='_')]]<<-unname(time[[k]])
    value
}
reference <- stamp('reference_verify_and_open',open_reference(spec$reference))
result <- stamp('processing',harmonise_file(spec$input,reference,spec$output,
                frequency_type='EAF',effect_type='beta',batch_size=as.integer(spec$batch_rows)))
if(!identical(result$implementation_version,spec$implementation))stop('Wrong normalization implementation')
if(result$source_rows!=spec$rows||result$normalized_rows!=spec$rows)stop('Row count mismatch')
metrics$source_rows<-result$source_rows
metrics$normalized_rows<-result$normalized_rows
metrics$status_counts<-result$status_counts
metrics$implementation<-result$implementation_version
metrics$reference_format<-result$core_version
metrics$reference_version<-result$reference_version
metrics$finished_utc<-format(Sys.time(),tz='UTC',usetz=TRUE)
jsonlite::write_json(metrics,spec$metrics,auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat('Completed',spec$version,result$source_rows,'source rows\n')
