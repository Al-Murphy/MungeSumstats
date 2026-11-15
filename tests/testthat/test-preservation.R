d <- tempfile('mss-'); dir.create(d)
writeLines(c('>1','CAAAAACGTACGT','>2','ACGT'), file.path(d,'tiny.fa'))
writeLines(c('chr\tpos\tref\talt\tid\taliases',
             '1\t8\tG\tA,T\trs8\told8;merged8',
             '1\t4\tAA\tA\trsDel\t.',
             '1\t10\tA\tC\trsDup\tshared',
             '1\t12\tG\tA\trsDup\tshared',
             '2\t2\tC\tG\tpalindrome\t.',
             '2\t1\tA\tC\tduplicate\t.',
             '1\t10\tA\tC\tduplicate\t.',
             '2\t3\tG\tA\tsplitThird\t.',
             '2\t3\tG\tC\tsplitThird2\t.'), file.path(d,'sites.tsv'))
test_that('immutable reference preparation records identity and rejects conflicts', {
  expect_error(prepare_reference(file.path(d,'tiny.fa'),file.path(d,'sites.tsv'),file.path(d,'ref')), 'build')
})
ref <- prepare_reference(file.path(d,'tiny.fa'),file.path(d,'sites.tsv'),file.path(d,'ref'),build='GRCh38',version='test-1',batch_size=2)
test_that('reference publication rejects input mutation during indexing', {
  mutable_fasta <- file.path(d,'mutable.fa');file.copy(file.path(d,'tiny.fa'),mutable_fasta)
  checked_prepare <- prepare_reference
  test_env <- new.env(parent=environment(prepare_reference))
  real_native <- get('mss_prepare_cpp',environment(prepare_reference))
  test_env$mss_prepare_cpp <- function(...) {
    real_native(...)
    cat('A',file=mutable_fasta,append=TRUE)
  }
  environment(checked_prepare) <- test_env
  expect_error(checked_prepare(mutable_fasta,file.path(d,'sites.tsv'),file.path(d,'mutated-ref'),build='GRCh38',version='changed'),'changed')
  expect_false(dir.exists(file.path(d,'mutated-ref')))
})
test_that('reference has immutable identity and source validation', {
  expect_identical(ref$build,'GRCh38')
  expect_identical(ref$version,'test-1')
  expect_identical(nchar(ref$checksums[['fasta']]),32L)
  expect_error(prepare_reference(file.path(d,'tiny.fa'),file.path(d,'sites.tsv'),file.path(d,'ref'),build='GRCh38',version='test-2'),'exists')
})
x <- data.frame(CHR=c('1','1','1','1','1','1','2','2',NA,'9','1','2'),
                BP=c('8','8','8','8','8','4','2','3',NA,'1','nonsense','1'),
                SNP=c('old8','rs8','rs8','rs8','rs8',NA,'palindrome','splitThird','duplicate','unknown','rs8',NA),
                A1=c('G','T','C','A','A','AA','C','A','A','N','G','A'),
                A2=c('T','G','A','C','T','A','G','G','C','<DEL>','T','C'),
                BETA=c('2','2','2','2','2','2','2','2','2',NA,'2','2'),
                FRQ=c('0','.2','.2','.2','.2','.2','.2','.2','.2',NA,'.2','0'),
                P=c('1e-1000',rep('0.5',11)),unknown=seq_len(12),check.names=FALSE)
y <- harmonise_sumstats(x,ref,frequency_type='EAF',batch_size=2)
test_that('source rows and raw numeric tokens survive all outcomes', {
  expect_identical(y$source,x)
  expect_identical(nrow(y$normalized),nrow(x))
  expect_identical(y$source$P[[1]],'1e-1000')
  expect_true(all(is.na(y$normalized$grch38_pos)))
  expect_identical(y$normalized$status,c('match','swap','complement','swapped_complement','alt_alt_contrast','ambiguous','ambiguous','swap','ambiguous','unmapped','invalid_position','match'))
})
test_that('effects, alleles and frequencies are conservative and per-allele', {
  expect_equal(y$normalized$effect[1:4],c(2,-2,2,-2))
  expect_equal(y$normalized$eaf[c(1,3,12)],c(0,.2,0))
  expect_true(all(is.na(y$normalized$eaf[c(2,4,7,8)])))
  expect_true(is.na(y$normalized$pos[[6]]))
  expect_true(is.na(y$normalized$ref[[6]]))
  expect_identical(y$normalized$alt[[1]],'T')
  expect_true(is.na(y$normalized$effect[[5]]))
  expect_equal(y$normalized$swapped[1:4],c(FALSE,TRUE,FALSE,TRUE))
  expect_equal(y$normalized$complemented[1:4],c(FALSE,FALSE,TRUE,TRUE))
})
test_that('split ALT contrasts remain distinct and generic opposite indel events stay ambiguous', {
  aa <- data.frame(CHR='2',BP=3,SNP=NA,A1='A',A2='C')
  expect_identical(harmonise_sumstats(aa,ref)$normalized$status,'alt_alt_contrast')
  shifted <- data.frame(CHR='1',BP=5,SNP=NA,A1='AA',A2='A')
  expect_identical(harmonise_sumstats(shifted,ref)$normalized$status,'ambiguous')
  expect_true(is.na(harmonise_sumstats(shifted,ref)$normalized$pos))
})
test_that('duplicate and unsupported reference rows retain third-allele context', {
  duplicate_sites <- file.path(d,'duplicate-sites.tsv')
  writeLines(c('chr\tpos\tref\talt\tid\taliases',
               '1\t8\tG\tA\tbiallelic\tshared-id',
               '1\t8\tG\tA,<DEL>\tmultiallelic\tshared-id',
               '1\t10\tA\tC\tsnp\t.',
               '1\t10\tA\t<DUP>\tsv\t.'),duplicate_sites)
  duplicated_ref <- prepare_reference(file.path(d,'tiny.fa'),duplicate_sites,file.path(d,'duplicate-ref'),build='GRCh38',version='duplicate-context')
  contrast <- data.frame(CHR='1',BP=8,SNP='multiallelic',A1='A',A2='G',BETA=.3,FRQ=.6)
  normalized <- harmonise_sumstats(contrast,duplicated_ref,frequency_type='EAF')$normalized
  expect_identical(normalized$status,'swap')
  expect_identical(normalized$id,'multiallelic')
  expect_true(is.na(normalized$eaf))
  contrast$SNP <- NA_character_
  expect_true(is.na(harmonise_sumstats(contrast,duplicated_ref,frequency_type='EAF')$normalized$id))
  contrast$BP <- 10;contrast$SNP <- 'snp';contrast$A1 <- 'C';contrast$A2 <- 'A'
  expect_true(is.na(harmonise_sumstats(contrast,duplicated_ref,frequency_type='EAF')$normalized$eaf))
})
test_that('frequency declarations and effect scale are explicit', {
  expect_true(all(is.na(harmonise_sumstats(x,ref)$normalized$eaf)))
  expect_true(all(is.na(harmonise_sumstats(x,ref,frequency_type='MAF')$normalized$eaf)))
  expect_true(all(is.na(harmonise_sumstats(x,ref,frequency_type='unknown')$normalized$eaf)))
  expect_true(all(is.na(harmonise_sumstats(x,ref,effect_type='none')$normalized$effect)))
  expect_error(harmonise_sumstats(x,ref,effect_type='odds_ratio'),'arg')
})
test_that('batching, empty input and repeated harmonisation preserve alignment', {
  expect_identical(y$normalized,harmonise_sumstats(x,ref,frequency_type='EAF',batch_size=100)$normalized)
  expect_identical(y$normalized,harmonise_sumstats(y,ref,frequency_type='EAF')$normalized)
  expect_identical(y$source,harmonise_sumstats(y,ref,frequency_type='EAF')$source)
  expect_equal(nrow(harmonise_sumstats(x[FALSE,],ref)$normalized),0)
  expect_error(harmonise_sumstats(x,ref,batch_size=0),'positive')
  expect_error(harmonise_sumstats(x,ref,build='GRCh37'),'build')
})
test_that('a verified native reference handle can be reused across calls', {
  handle <- open_reference(ref)
  expect_s3_class(handle,'mss_reference_handle')
  expect_identical(harmonise_sumstats(x,handle,frequency_type='EAF')$normalized,y$normalized)
  expect_identical(harmonise_sumstats(x,handle,frequency_type='EAF',batch_size=1)$normalized,y$normalized)
})
test_that('file streaming matches resident batches and preserves associations', {
  input <- file.path(d,'input.tsv')
  write.table(x,input,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  out <- harmonise_file(input,ref,file.path(d,'stream'),batch_size=3)
  resident <- harmonise_sumstats(read.delim(input,colClasses='character',check.names=FALSE),ref)
  disk <- read.delim(file.path(out$path,'normalized.tsv'),check.names=FALSE)
  expect_identical(disk$status,resident$normalized$status)
  expect_equal(disk$source_row,seq_len(nrow(x)))
  expect_equal(disk$effect,resident$normalized$effect)
  expect_equal(out$source_rows,nrow(x))
  expect_true(all(disk$normalized_build=='GRCh38'))
  metadata <- read.delim(file.path(out$path,'metadata.tsv'),colClasses='character')
  expect_identical(metadata$value[metadata$key=='reference_build'],'GRCh38')
  expect_identical(metadata$value[metadata$key=='source_file'],'source.tsv')
  receipts <- read.delim(file.path(out$path,'checksums.tsv'),colClasses='character')
  expect_identical(unname(tools::md5sum(file.path(out$path,receipts$filename))),receipts$md5)
  expect_true(file.exists(file.path(out$path,'columns.tsv')))
  expect_true(file.exists(file.path(out$path,'reference_checksums.tsv')))
  expect_equal(sum(out$status_counts),nrow(x))
  expect_identical(unname(tools::md5sum(input)),unname(tools::md5sum(file.path(out$path,'source-artifact'))))
  expect_identical(readLines(input),readLines(file.path(out$path,'source.tsv')))
  gz <- file.path(d,'input.tsv.gz'); con <- gzfile(gz,'wt');writeLines(readLines(input),con);close(con)
  outgz <- harmonise_file(gz,ref,file.path(d,'stream-gz'),batch_size=2)
  expect_identical(readLines(file.path(out$path,'normalized.tsv')),readLines(file.path(outgz$path,'normalized.tsv')))
  expect_error(harmonise_file(input,ref,out$path),'exists')
})
test_that('malformed or truncated input cannot publish an edition', {
  bad <- file.path(d,'bad.tsv')
  writeLines(c('CHR\tBP\tSNP\tA1\tA2','1\t8\trs8\tG\tT','1\t8\trs8'),bad)
  expect_error(harmonise_file(bad,ref,file.path(d,'failed'),batch_size=1),'fields')
  expect_false(dir.exists(file.path(d,'failed')))
  con <- file(bad,'wb');writeChar('CHR\tBP\n1\t8',con,eos=NULL);close(con)
  expect_error(harmonise_file(bad,ref,file.path(d,'truncated'),batch_size=1),'read|line')
  expect_false(dir.exists(file.path(d,'truncated')))
  good_gz <- file.path(d,'complete.tsv.gz'); con <- gzfile(good_gz,'wt');writeLines(c('CHR\tBP','1\t8'),con);close(con)
  bytes <- readBin(good_gz,'raw',n=file.info(good_gz)$size)
  bad_gz <- file.path(d,'broken.tsv.gz');writeBin(bytes[seq_len(length(bytes)-8)],bad_gz)
  expect_error(harmonise_file(bad_gz,ref,file.path(d,'broken-gzip'),batch_size=1),'gzip|read|compressed')
  expect_false(dir.exists(file.path(d,'broken-gzip')))
})
test_that('reference payload mutation fails verification', {
  cat('\n',file=file.path(d,'ref','records'),append=TRUE)
  expect_error(harmonise_sumstats(x,ref),'checksum')
})
unlink(d,recursive=TRUE)
