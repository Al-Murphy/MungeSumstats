root <- '/private/tmp/gwas-harvest-stage2/review-mss-context-format3-stage7'
base <- '/private/tmp/gwas-harvest-stage2'
lib <- file.path(base,'results/mss-context-format3-stage7/installed-lib')
.libPaths(c(lib,file.path(base,'Rlib'),'/private/tmp/harvest-Rlib',.libPaths()))
suppressPackageStartupMessages(library(MungeSumstats));library(testthat)
stopifnot(normalizePath(find.package('MungeSumstats'))==normalizePath(file.path(lib,'MungeSumstats')))
ns <- asNamespace('MungeSumstats')
stopifnot(getNativeSymbolInfo('_MungeSumstats_mss_prepare_cpp',PACKAGE='MungeSumstats')$numParameters==5L)
e <- new.env(parent=ns)
results <- test_dir(file.path(root,'source/tests/testthat'),filter='reference-format3',env=e,reporter='summary',stop_on_failure=TRUE)
x <- as.data.frame(results)
stopifnot(nrow(x)==7L,sum(x$nb)==65L,sum(x$failed)==0L,sum(x$error)==0L)
dir <- file.path(root,'installed-controls');dir.create(dir)
setwd(dir)
checks <- character()
check <- function(value,label){stopifnot(isTRUE(value));checks<<-c(checks,label)}
refuse <- function(expr,label){result<-tryCatch({force(expr);FALSE},error=function(e)TRUE);check(result,label)}
writeLines(c('>1','CAAAAACGTACGT'), 'reference.fa')
sites <- c('chr\tpos\tref\talt\tid\taliases','1\t8\tG\tA,T\tr\told;同義;é;é','1\t4\tAA\tA\tdel\told_del','1\t10\tN\t<DEL>,*\tunsupported\tunsupported_alias')
writeLines(sites,'sites.tsv',useBytes=TRUE)
native <- get('mss_prepare_cpp',ns)
native('reference.fa','sites.tsv','native-four',2L)
check(identical(readLines('native-four/core-version'),'mss-preservation-2'),'native R wrapper old four arguments defaults to2')
refs <- list()
for (build in c('GRCh37','GRCh38'))for(format in 2:3) {
 key<-paste(build,format)
 refs[[key]]<-prepare_reference('reference.fa','sites.tsv',paste0('ref-',build,'-',format),build=build,version='typed-test',reference_format=paste0('mss-preservation-',format),contig_aliases=c(chr1='1'))
}
default <- prepare_reference('reference.fa','sites.tsv','default',build='GRCh37',version='typed-test',contig_aliases=c(chr1='1'))
for (file in c('core-version','fasta.fai','keys','offsets','records','rejected.tsv','source.tsv','contig_aliases.tsv'))
 check(identical(readBin(file.path(default$path,file),'raw',file.info(file.path(default$path,file))$size),readBin(file.path(refs[['GRCh37 2']]$path,file),'raw',file.info(file.path(refs[['GRCh37 2']]$path,file))$size)),paste('default2 explicit2 exact',file))
check(identical(readBin('sites.tsv','raw',file.info('sites.tsv')$size),readBin(file.path(refs[['GRCh37 3']]$path,'source.tsv'),'raw',file.info('sites.tsv')$size)),'complete unsupported IDs and source bytes retained')
q<-data.frame(CHR='chr1',BP=c(8,8,4),SNP=c('同義','old','old_del'),A1=c('G','T','AA'),A2=c('T','G','A'),BETA=c('.3','-.3','.125'),FRQ=c('.2','.8','.25'),extra='0001',stringsAsFactors=FALSE)
y<-lapply(2:3,function(format)harmonise_sumstats(q,refs[[paste('GRCh37',format)]],frequency_type='EAF'))
check(identical(y[[1]]$normalized,y[[2]]$normalized),'public full normalized2/3 exact with contig and historical aliases')
for(i in 1:2)check(identical(y[[i]]$source,q),paste('source all lexical fields exact',i))
check(identical(y[[2]]$provenance$implementation_version,'mss-preservation-7'),'science version7 unchanged')
check(identical(y[[2]]$provenance$core_version,'mss-preservation-3'),'format3 data identity receipt')
check(identical(y[[2]]$provenance$reference_context_order,'semantic-s-context-v1'),'format3 semantic context order receipt')
writeLines(c('chain 100 1 13 + 0 13 1 13 + 0 13 1','13'),'identity.chain')
old <- new.env(parent=ns);sys.source(file.path(base,'review-mss-preservation-v7-stage7/source/R/projection.R'),old)
expectation<-NULL
for(sf in 2:3)for(tf in 2:3) {
 p<-prepare_projection('identity.chain',refs[[paste('GRCh37',sf)]],refs[[paste('GRCh38',tf)]],'GRCh37','GRCh38')
 projected<-harmonise_sumstats(q,refs[[paste('GRCh37',sf)]],frequency_type='EAF',projection=p)$normalized
 if(is.null(expectation))expectation<-projected else check(identical(projected,expectation),paste('projection full scientific fields exact',sf,tf))
 check(identical(p$receipt$reference_format,if(sf==tf)paste0('mss-preservation-',sf) else 'mixed'),paste('projection common format',sf,tf))
 check(identical(p$receipt$normalization_implementation,'mss-preservation-7'),paste('projection normalizer identity',sf,tf))
 if(sf==2&&tf==2) {
  previous<-old$prepare_projection('identity.chain',refs[['GRCh37 2']],refs[['GRCh38 2']],'GRCh37','GRCh38')
  check(identical(p$receipt,previous$receipt),'2to2 complete projection receipt unchanged from SCI7')
 } else {
  check(identical(p$receipt$source_reference_format,paste0('mss-preservation-',sf)),paste('projection source format',sf,tf))
  check(identical(p$receipt$target_reference_format,paste0('mss-preservation-',tf)),paste('projection target format',sf,tf))
 }
}
clone <- function(name){dir.create(name);stopifnot(all(file.copy(list.files(refs[['GRCh37 3']]$path,full.names=TRUE),name)));normalizePath(name)}
for(field in c('record_layout','reference_context_order','core_version')) {
 dest<-clone(paste0('bad-',field));m<-readRDS(file.path(dest,'manifest.rds'));m[[field]]<-'wrong';saveRDS(m,file.path(dest,'manifest.rds'))
 refuse(open_reference(dest),paste('public manifest mismatch refuses',field))
}
dest<-clone('bad-marker');writeLines('mss-preservation-2',file.path(dest,'core-version'));refuse(open_reference(dest),'public marker mismatch refuses')
dest<-clone('bad-payload');b<-readBin(file.path(dest,'records'),'raw',file.info(file.path(dest,'records'))$size);b[[length(b)-1L]]<-as.raw(88L);writeBin(b,file.path(dest,'records'));refuse(open_reference(dest),'public checksum refuses corrupt unqueried payload')
b<-readBin('sites.tsv','raw',file.info('sites.tsv')$size);b[[length(b)-2L]]<-as.raw(255L);writeBin(b,'invalid.tsv')
refuse(prepare_reference('reference.fa','invalid.tsv','invalid-published',build='test',version='v',reference_format='mss-preservation-3'),'public invalid UTF8 source refuses')
check(!file.exists('invalid-published'),'invalid source publishes no output')
write.table(data.frame(check=checks,passed=TRUE),file.path(root,'installed-controls.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(sprintf('installed API suite: %d blocks / %d assertions; independent installed controls: %d',nrow(x),sum(x$nb),length(checks)),file.path(root,'installed-summary.txt'))
cat(readLines(file.path(root,'installed-summary.txt')),'\n')
