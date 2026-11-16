root <- normalizePath(commandArgs(trailingOnly=TRUE)[[1]])
.libPaths(c('/private/tmp/harvest-Rlib', .libPaths()))
Rcpp::sourceCpp(file.path(root,'source/src/preservation.cpp'))
source(file.path(root,'source/R/preservation.R'))
source(file.path(root,'source/R/projection.R'))
out <- file.path(root,'public-R');dir.create(out)
checks <- 0L
check <- function(x) {checks <<- checks+1L;stopifnot(x)}
fasta <- file.path(out,'tiny.fa');writeLines(c('>X','GACAAAAA'),fasta)
make <- function(name,row,build='GRCh37') {
  site <- file.path(out,paste0(name,'.tsv'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases',row),site)
  prepare_reference(fasta,site,file.path(out,name),build,'independent-stage6')
}
plain <- make('plain','X\t2\tA\tC,G,T\tr\told_r')
padded <- make('padded','X\t1\tGA\tGC,GG,GT\tr\told_r')
input <- data.frame(CHR='X',BP=2,SNP='old_r',A1=c('A','G','T','C'),
                    A2=c('G','A','C','T'),BETA='-0.625',FRQ='0.137',
                    P='1e-999',token=c('00001','NA','source-token','1.000'),stringsAsFactors=FALSE)
a <- harmonise_sumstats(input,plain,frequency_type='EAF')
b <- harmonise_sumstats(input,padded,frequency_type='EAF')
check(identical(a$source,input));check(identical(b$source,input))
check(identical(a$normalized$status,rep('ambiguous',4)))
check(all(is.na(a$normalized[c('chr','pos','ref','alt','id','effect','eaf','swapped','complemented')])))
# The reanchored correction must now agree with the biological minimal-site decision.
check(identical(b$normalized$status,rep('ambiguous',4)))
check(all(is.na(b$normalized[c('chr','pos','ref','alt','id','effect','eaf','swapped','complemented')])))
for (kind in c('MAF','unknown')) {
  z <- harmonise_sumstats(input,padded,frequency_type=kind)
  check(identical(z$normalized$status,b$normalized$status))
  check(all(is.na(z$normalized$eaf)))
}
background <- make('background','X\t2\tAC\tGC,GT,AT\ts\told_s')
background_input <- input
background_input$BP <- 3; background_input$SNP <- 'old_s'
background_input$A1 <- c('C','T','G','A'); background_input$A2 <- c('T','C','A','G')
for (kind in c('EAF','MAF','unknown')) {
  z <- harmonise_sumstats(background_input,background,frequency_type=kind)
  check(identical(z$source,background_input))
  check(identical(z$normalized$status,rep('ambiguous',4)))
  check(all(is.na(z$normalized[c('chr','pos','ref','alt','id','effect','eaf','swapped','complemented')])))
}
# Explicit measured source REF/ALT and strand legitimately select one hypothesis.
columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',
                source_ref='REF',source_alt='ALT',source_strand='STRAND')
declared <- input
declared$REF <- c('A','A','T','T');declared$ALT <- c('G','G','C','C');declared$STRAND <- c('+','+','-','-')
trusted <- harmonise_sumstats(declared,padded,columns=columns,frequency_type='EAF')
check(identical(trusted$source,declared))
check(identical(trusted$normalized$status,c('match','swap','complement','swapped_complement')))
check(identical(trusted$normalized$effect,c(-.625,.625,-.625,.625)))
declared$STRAND <- 'unknown'
check(identical(harmonise_sumstats(declared,padded,columns=columns)$normalized$status,rep('source_strand_unknown',4)))
check(identical(a$provenance$implementation_version,'mss-preservation-7'))
check(identical(a$provenance$core_version,'mss-preservation-2'))
target <- make('target','X\t2\tA\tC,G,T\tr\told_r','GRCh38')
chain <- file.path(out,'identity.chain')
writeLines(c('chain 100 X 8 + 0 8 X 8 + 0 8 9','8'),chain)
projection <- prepare_projection(chain,plain,target,'GRCh37','GRCh38')
check(identical(projection$receipt$normalization_implementation,'mss-preservation-7'))
check(identical(projection$receipt$implementation_version,'mss-projection-2'))
check(identical(projection$receipt$reference_format,'mss-preservation-2'))
input_path <- file.path(out,'input.tsv');write.table(input,input_path,sep='\t',row.names=FALSE,quote=TRUE)
edition <- harmonise_file(input_path,plain,file.path(out,'edition'),frequency_type='EAF',batch_size=1)
check(identical(unname(tools::md5sum(input_path)),unname(tools::md5sum(file.path(edition$path,'source-artifact')))))
check(edition$source_rows==4 && edition$normalized_rows==4)
manifest <- readRDS(file.path(edition$path,'manifest.rds'))
check(identical(manifest$implementation_version,'mss-preservation-7'))
check(identical(manifest$core_version,'mss-preservation-2'))
metadata <- read.delim(file.path(edition$path,'metadata.tsv'),colClasses='character')
check(identical(metadata$value[metadata$key=='implementation_version'],'mss-preservation-7'))
write.table(data.frame(result='PASS: padded and shared-background ambiguity, contracts, receipts, streaming verified',assertions=checks),
            file.path(root,'public-R-results.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write.table(a$normalized,file.path(out,'plain-public-output.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write.table(b$normalized,file.path(out,'padded-public-output.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
cat('Independent public R assertions:',checks,'PASS; corrected padded decisions and shared-background ambiguity verified\n')

context_fixture <- function(name, sequence, rows) {
  dir <- file.path(out,name);dir.create(dir)
  fa <- file.path(dir,'tiny.fa');writeLines(c('>Z',sequence),fa)
  site <- file.path(dir,'sites.tsv');writeLines(c('chr\tpos\tref\talt\tid\taliases',rows),site)
  prepare_reference(fa,site,file.path(dir,'index'),'GRCh37',name)
}
empty <- function(x,status) {
  check(identical(x$status,rep(status,nrow(x))))
  check(all(is.na(x[c('chr','pos','ref','alt','id','effect','eaf','swapped','complemented')])) )
}
make_query <- function(a,b,pos,id='old_r') data.frame(CHR='Z',BP=pos,SNP=id,A1=a,A2=b,
    BETA='0.625',FRQ='0.137',P='1e-1000',token='00001',stringsAsFactors=FALSE)
unique <- context_fixture('unknown-anchor','CCAGG','Z\t1\tCCAG\tCCTG\tu\told_u')
q <- make_query(c('CCAG','CCTG','CTGG','CAGG'),c('CCTG','CCAG','CAGG','CTGG'),NA_real_,'old_u')
for(kind in c('EAF','MAF','unknown')) {
  z <- harmonise_sumstats(q,unique,frequency_type=kind)
  check(identical(z$source,q));empty(z$normalized,'source_coordinates_missing')
  check(identical(z$provenance$implementation_version,'mss-preservation-7'))
}
known <- make_query('CCAG','CCTG',c(1,2),'old_u')
z <- harmonise_sumstats(known,unique,frequency_type='EAF')
check(identical(z$normalized$status,c('match','swapped_complement')))
check(identical(z$normalized$effect,c(.625,-.625)))
check(isTRUE(all.equal(z$normalized$eaf,c(.137,1-.137))))
for(spec in list(
    list(name='repeat-context',seq='CAAAAAAA',row='Z\t1\tCAA\tGA,GAAA,CAAAA\tr\told_r',
         a=c('A','AAA','T','TTT'),b=c('AAA','A','TTT','T'),pos=4),
    list(name='inserted-background',seq='AAAAAA',row='Z\t1\tAA\tGCA,GCG,AG\tr\told_r',
         a=c('A','G','T','C'),b=c('G','A','C','T'),pos=2),
    list(name='unknown-sequence',seq='AAAAAA',row='Z\t1\tA\tC,<DEL>\tr\told_r',
         a='A',b='C',pos=1))) {
  reference <- context_fixture(spec$name,spec$seq,spec$row)
  source <- make_query(spec$a,spec$b,spec$pos)
  for(kind in c('EAF','MAF','unknown')) {
    z <- harmonise_sumstats(source,reference,frequency_type=kind)
    check(identical(z$source,source));empty(z$normalized,'source_allele_context_required')
  }
  if(spec$name=='repeat-context') {
    trusted <- make_query('A','AAA',4)
    trusted$REF <- 'A';trusted$ALT <- 'AAA';trusted$STRAND <- '+'
    z <- harmonise_sumstats(trusted,reference,columns=columns,frequency_type='EAF')
    check(identical(z$source,trusted));check(identical(z$normalized$status,'match'))
    check(identical(z$normalized$effect,.625));check(identical(z$normalized$eaf,.137))
  }
  path <- file.path(out,paste0(spec$name,'.tsv'));write.table(source,path,sep='\t',quote=TRUE,row.names=FALSE)
  for(batch in c(1,3)) {
    edition <- harmonise_file(path,reference,file.path(out,paste0(spec$name,'-edition-',batch)),frequency_type='EAF',batch_size=batch)
    check(edition$source_rows==nrow(source) && edition$normalized_rows==nrow(source))
    check(identical(unname(tools::md5sum(path)),unname(tools::md5sum(file.path(edition$path,'source-artifact')))))
    normalized <- read.delim(file.path(edition$path,'normalized.tsv'),check.names=FALSE)
    empty(normalized,'source_allele_context_required')
    check(identical(readRDS(file.path(edition$path,'manifest.rds'))$implementation_version,'mss-preservation-7'))
  }
}
missing_path <- file.path(out,'missing-position.tsv');write.table(q,missing_path,sep='\t',quote=TRUE,row.names=FALSE)
edition <- harmonise_file(missing_path,unique,file.path(out,'missing-position-edition'),frequency_type='EAF',batch_size=1)
check(edition$source_rows==4 && edition$normalized_rows==4)
empty(read.delim(file.path(edition$path,'normalized.tsv'),check.names=FALSE),'source_coordinates_missing')
check(identical(unname(tools::md5sum(missing_path)),unname(tools::md5sum(file.path(edition$path,'source-artifact')))))
budget <- context_fixture('budget','AAAAA',paste0('Z\t1\tA\tG,',paste(rep('N',70),collapse=','),'\tr\t.'))
input <- make_query('A','G',1,'r');path <- file.path(out,'budget-input.tsv');write.table(input,path,sep='\t',row.names=FALSE,quote=TRUE)
destination <- file.path(out,'budget-refused-edition')
failed <- tryCatch({harmonise_file(path,budget,destination,candidate_limit=2);FALSE},
                   error=function(e)grepl('context closure limit exceeded',conditionMessage(e)))
check(failed);check(!file.exists(destination));check(file.exists(path));check(file.exists(file.path(budget$path,'records')))
write.table(data.frame(result='PASS: independent v7 API/source/stream/receipt/fail-closed checks',assertions=checks),
            file.path(root,'public-R-results.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
cat('Independent final public R assertions:',checks,'PASS\n')
