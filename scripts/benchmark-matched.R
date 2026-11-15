#!/usr/bin/env Rscript
# Matched QC benchmark; see scripts/benchmark-matched.md for its explicit scope.
a <- commandArgs(TRUE)
if (length(a) != 7L) stop('mode input reference output metrics cache threads')
names(a) <- c('mode','input','reference','output','metrics','cache','threads')
stage_lib <- Sys.getenv('MSS_STAGE_RLIB','/private/tmp/gwas-harvest-stage2/Rlib')
.libPaths(c(if(a[['mode']]=='legacy') file.path(stage_lib,'legacy') else Sys.getenv('MSS_NEW_RLIB',stage_lib), stage_lib,
            '/private/tmp/harvest-Rlib',.libPaths()))
library(MungeSumstats)
library(data.table)
setDTthreads(as.integer(a[['threads']]))
metrics <- list(mode=a[['mode']], cache=a[['cache']], input=normalizePath(a[['input']]),
                reference=normalizePath(a[['reference']]), R_version=R.version.string,
                threads=a[['threads']], package_version=as.character(packageVersion('MungeSumstats')),
                os_cache='uncontrolled; reference-cold is application state, not flushed OS cache',
                legacy_commit='b4a2d21', new_initial_commit='d57c8f9')
metrics$loaded_dll_md5 <- if(a[['mode']]=='new')
    unname(tools::md5sum(getLoadedDLLs()[['MungeSumstats']][['path']])) else 'legacy_has_no_native_DLL'
stamp <- function(name, expr) {
    time <- system.time(value <- force(expr))
    for (part in c('elapsed','user.self','sys.self'))
        metrics[[paste(name,part,sep='_')]] <<- unname(time[[part]])
    value
}
if (a[['mode']]=='new') {
    ref <- stamp('reference_open',open_reference(a[['reference']]))
    if (a[['cache']]=='warm') {
        # Prime the complete identical workload, preserving the verified handle.
        invisible(stamp('prime',harmonise_file(a[['input']],ref,paste0(a[['output']],'.prime'),
                              frequency_type='EAF',batch_size=100000)))
        metrics$prime_scope <- 'complete identical input, same verified handle'
        gc()
    }
    result <- stamp('processing',harmonise_file(a[['input']],ref,a[['output']],
                       frequency_type='EAF',batch_size=100000))
    metrics$source_rows <- result$source_rows
    metrics$retained_rows <- result$normalized_rows
    metrics$status_counts <- paste(names(result$status_counts),result$status_counts,sep=':',collapse=';')
    metrics$reference_version <- result$reference_version
    metrics$core_version <- result$core_version
    metrics$implementation_version <- if(is.null(result$implementation_version)) 'mss-preservation-2' else result$implementation_version
} else if (a[['mode']]=='legacy') {
    # This adapter replaces only the legacy reference acquisition, keeping all
    # format_sumstats QC functions from the installed b4a2d21 package intact.
    load_panel <- function() {
        sites <- fread(file.path(a[['reference']],'source.tsv'))
        map <- Biostrings::IUPAC_CODE_MAP
        # All records in the matched panel are biallelic SNPs. Construct the
        # IUPAC mapping once, not once per reference row.
        pairs <- vapply(map,function(s) paste(sort(strsplit(s,'')[[1]]),collapse=''),character(1))
        codes <- setNames(names(map),pairs)
        pair <- ifelse(sites$ref < sites$alt,paste0(sites$ref,sites$alt),paste0(sites$alt,sites$ref))
        ambiguity <- unname(codes[pair])
        if(anyNA(ambiguity)) stop('Matched reference must contain biallelic SNPs')
        out <- data.table(SNP=sites$id,seqnames=sites$chr,pos=sites$pos,
                          alleles_as_ambig=ambiguity,ref_allele=sites$ref,alt_alleles=sites$alt)
        setkey(out,SNP)
        out
    }
    panel <- if(a[['cache']]=='warm') stamp('reference_open',load_panel()) else NULL
    ns <- asNamespace('MungeSumstats')
    adapter <- function(snps,...) {
        value <- if(is.null(panel)) stamp('reference_open',load_panel()) else panel
        out <- value[unique(snps),nomatch=NULL]
        setkey(out,SNP)
        out
    }
    unlockBinding('load_ref_genome_data',ns)
    assign('load_ref_genome_data',adapter,envir=ns)
    lockBinding('load_ref_genome_data',ns)
    metrics$reference_version <- readRDS(file.path(a[['reference']],'manifest.rds'))$version
    legacy_call <- function(output) format_sumstats(a[['input']],ref_genome='GRCh37',dbSNP=144,
                  on_ref_genome=FALSE,infer_eff_direction=FALSE,check_dups=FALSE,
                  strand_ambig_filter=FALSE,frq_is_maf=FALSE,sort_coordinates=FALSE,
                  nThread=as.integer(a[['threads']]),save_path=output,return_data=TRUE)
    if(a[['cache']]=='warm') {
        prime <- stamp('prime',legacy_call(paste0(a[['output']],'.prime.tsv')))
        rm(prime);gc()
        metrics$prime_scope <- 'complete identical input, same loaded reference panel'
    }
    result <- stamp('processing',legacy_call(a[['output']]))
    metrics$retained_rows <- nrow(result)
    metrics$implementation_version <- 'legacy-b4a2d21'
    expected_rows <- Sys.getenv('MSS_EXPECTED_ROWS')
    metrics$source_rows <- if(nzchar(expected_rows)) as.numeric(expected_rows) else
        sum(count.fields(a[['input']],sep='\t',quote='')>0)-1
    metrics$source_count_origin <- if(nzchar(expected_rows)) 'acquisition/generator receipt; full file MD5 recorded' else 'count.fields full scan'
    metrics$status_counts <- 'legacy_has_no_row_statuses'
    metrics$adapter <- 'reference loader only; all legacy QC unchanged; on_ref_genome disabled for all-known matched SNPs; duplicates allowed'
} else stop('mode must be new or legacy')
metrics$input_md5 <- unname(tools::md5sum(a[['input']]))
metrics$created_utc <- format(Sys.time(),tz='UTC',usetz=TRUE)
write.table(data.frame(key=names(metrics),value=vapply(metrics,as.character,character(1))),
            a[['metrics']],sep='\t',quote=TRUE,row.names=FALSE)
cat('Completed',a[['mode']],metrics$source_rows,metrics$retained_rows,'rows\n')
