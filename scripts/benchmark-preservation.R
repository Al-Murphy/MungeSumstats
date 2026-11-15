#!/usr/bin/env Rscript
# Caller-supplied asset benchmark. No performance or RSS conclusion is implied.
usage <- function() cat(paste(
  'Usage: Rscript scripts/benchmark-preservation.R',
  '  --input FILE --reference DIRECTORY --output NEW_DIRECTORY',
  '  --metrics NEW_TSV --rows EXPECTED_SOURCE_ROWS',
  '  [--batch-size 100000] [--frequency-type unknown|EAF|MAF]',
  '  [--effect-type beta|log_odds|none] [--build LABEL]',
  '  [--columns chr=CHR,pos=BP,id=SNP,other=A1,effect_allele=A2,effect=BETA,eaf=FRQ]',
  '  [--delimiter tab|comma] [--cache-state caller-declared-description]',
  '  [--source-root MSS_CHECKOUT]',
  '',
  'Reference opening/hash verification and file harmonisation are timed separately.',
  'The timed file path includes snapshotting, parsing, normalization, output and hashes.',
  'Cache state is only caller-declared; this script does not flush OS/reference caches.',
  'R allocations are not RSS. To collect RSS externally, wrap the command with',
  '  Linux: /usr/bin/time -v Rscript ...',
  '  macOS: /usr/bin/time -l Rscript ...',
  'Use an OS process-memory monitor on Windows; save its output beside metrics.',
  'No RSS, speedup, cold/warm equivalence or production target is inferred here.',
  sep='\n'), '\n')
argv <- commandArgs(trailingOnly=TRUE)
if (!length(argv) || '--help' %in% argv) { usage();quit(status=0L) }
allowed <- c('input','reference','output','metrics','rows','batch-size',
             'frequency-type','effect-type','build','columns','delimiter',
             'cache-state','source-root')
if (length(argv) %% 2L) stop('Arguments must be --name value pairs; use --help')
options <- list()
for (i in seq.int(1L,length(argv),by=2L)) {
  flag <- argv[[i]]
  if (!startsWith(flag,'--')) stop('Expected --name: ',flag)
  name <- substring(flag,3L)
  if (!name %in% allowed || !is.null(options[[name]])) stop('Unknown or duplicate argument: ',flag)
  options[[name]] <- argv[[i+1L]]
}
for (name in c('input','reference','output','metrics','rows'))
  if (is.null(options[[name]]) || !nzchar(options[[name]])) stop('Missing --',name)
expected <- suppressWarnings(as.numeric(options$rows))
if (length(expected)!=1L || !is.finite(expected) || expected<0 || expected!=floor(expected) || expected>2^53-1)
  stop('--rows must be an exact non-negative expected row count')
if (file.exists(options$metrics)) stop('--metrics path already exists')
if (file.exists(options$output)) stop('--output path already exists')
metrics_path <- file.path(normalizePath(dirname(options$metrics),mustWork=TRUE),basename(options$metrics))
get_option <- function(name,default) if(is.null(options[[name]])) default else options[[name]]
batch_size <- suppressWarnings(as.numeric(get_option('batch-size','100000')))
frequency_type <- match.arg(get_option('frequency-type','unknown'),c('unknown','EAF','MAF'))
effect_type <- match.arg(get_option('effect-type','beta'),c('beta','log_odds','none'))
delimiter_name <- match.arg(get_option('delimiter','tab'),c('tab','comma'))
cache_state <- get_option('cache-state','unspecified-uncontrolled')
columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ')
if (!is.null(options$columns)) {
  entries <- strsplit(options$columns,',',fixed=TRUE)[[1L]]
  fields <- strsplit(entries,'=',fixed=TRUE)
  if (any(lengths(fields)!=2L)) stop('--columns requires field=source_column pairs')
  keys <- vapply(fields,`[[`,character(1),1L)
  if (anyDuplicated(keys) || any(!keys %in% names(columns))) stop('Unknown/duplicate --columns field')
  columns <- stats::setNames(lapply(fields,`[[`,2L),keys)
}
test_lib <- Sys.getenv('MSS_TEST_RLIB')
if (nzchar(test_lib)) .libPaths(c(test_lib,.libPaths()))
if (!is.null(options[['source-root']])) {
  source_root <- normalizePath(options[['source-root']],mustWork=TRUE)
  Rcpp::sourceCpp(file.path(source_root,'src/preservation.cpp'))
  source(file.path(source_root,'R/preservation.R'))
} else {
  if (!requireNamespace('MungeSumstats',quietly=TRUE))
    stop('Install MungeSumstats or provide --source-root for standalone Rcpp compilation')
  open_reference <- getExportedValue('MungeSumstats','open_reference')
  harmonise_file <- getExportedValue('MungeSumstats','harmonise_file')
}
reference_time <- system.time(handle <- open_reference(options$reference))
processing_time <- system.time(result <- harmonise_file(
  options$input,handle,options$output,columns=columns,
  build=options$build,frequency_type=frequency_type,effect_type=effect_type,
  delimiter=if(delimiter_name=='tab') '\t' else ',',batch_size=batch_size))
count_valid <- identical(as.double(result$source_rows),as.double(expected)) &&
               identical(as.double(result$normalized_rows),as.double(expected))
metadata <- c(schema_version='1',benchmark='mss-file-preservation',
              expected_rows=sprintf('%.17g',expected),source_rows=sprintf('%.17g',result$source_rows),
              normalized_rows=sprintf('%.17g',result$normalized_rows),
              count_validation=if(count_valid) 'pass' else 'fail',
              batch_size=sprintf('%.17g',batch_size),frequency_type=frequency_type,
              effect_type=effect_type,delimiter=delimiter_name,
              reference_open_elapsed_seconds=sprintf('%.17g',reference_time[['elapsed']]),
              harmonise_elapsed_seconds=sprintf('%.17g',processing_time[['elapsed']]),
              harmonise_user_seconds=sprintf('%.17g',processing_time[['user.self']]),
              harmonise_system_seconds=sprintf('%.17g',processing_time[['sys.self']]),
              rss_measured='false',cache_state=cache_state,
              cache_caveat='Caller declared only; OS/reference caches were not controlled or flushed.',
              timing_scope='Reference opening separate; file timing includes snapshot, parse, native work, writes and hashes.',
              core_version=result$core_version,reference_build=result$reference_build,
              reference_version=result$reference_version,reference_path=result$reference_path,
              input_path=normalizePath(options$input,mustWork=TRUE),edition_path=result$path,
              input_checksum=result$input_checksum,checksum_algorithm=result$checksum_algorithm,
              projection=result$projection,R_version=R.version.string,
              platform=R.version$platform,created_utc=format(Sys.time(),tz='UTC',usetz=TRUE))
metadata <- c(metadata,stats::setNames(unname(result$checksums),paste0('output_md5:',names(result$checksums))),
              stats::setNames(unname(result$reference_checksums),paste0('reference_md5:',names(result$reference_checksums))))
utils::write.table(data.frame(key=names(metadata),value=unname(metadata)),metrics_path,
                   sep='\t',quote=TRUE,row.names=FALSE,qmethod='double')
cat('Metrics:',metrics_path,'\nEdition:',result$path,'\n')
cat('Source/output/expected rows:',result$source_rows,result$normalized_rows,expected,'\n')
cat('Elapsed file seconds:',processing_time[['elapsed']],'; RSS not measured; cache state uncontrolled.\n')
if (!count_valid) stop('Audited source/output row counts differ from --rows; metrics record count_validation=fail')
