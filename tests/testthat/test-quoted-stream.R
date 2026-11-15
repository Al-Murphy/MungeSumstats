test_that('quoted CSV preserves bytes and decoded token spelling across bounded chunks', {
  d <- tempfile(); dir.create(d)
  writeLines(c('>1','ACGT'),file.path(d,'f.fa'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t2\tC\tA\trs2\t.'),file.path(d,'v.tsv'))
  ref <- prepare_reference(file.path(d,'f.fa'),file.path(d,'v.tsv'),file.path(d,'ref'),'GRCh37','toy')
  bytes <- charToRaw('CHR,BP,SNP,A1,A2,BETA,FRQ,NOTE,EMPTY\r\n1,2,rs2,C,A,0.12345678901234566,.20,"one, two\r\nsaid ""hello""",\r\n1,2,rs2,C,A,1e-400,.3,tail,')
  input <- file.path(d,'i.csv'); writeBin(bytes,input)
  out <- harmonise_file(input,ref,file.path(d,'out'),delimiter=',',batch_size=1,batch_bytes=32)
  expect_equal(out$source_rows,2)
  expect_identical(readBin(file.path(out$path,'source-artifact'),'raw',length(bytes)),bytes)
  expect_identical(readBin(file.path(out$path,'source.csv'),'raw',length(bytes)),bytes)
  parsed <- suppressWarnings(read.csv(file.path(out$path,'source.csv'),colClasses='character',na.strings=NULL))
  expect_identical(parsed$BETA,c('0.12345678901234566','1e-400'))
  expect_identical(parsed$NOTE,c('one, two\nsaid "hello"','tail'))
  norm <- read.delim(file.path(out$path,'normalized.tsv'))
  expect_equal(norm$source_row,1:2)
  expect_equal(norm$effect[1],0.12345678901234566,tolerance=0)
  for (bad in c('CHR,BP\n1,"2','CHR,BP\n1,2"','CHR,BP\n1,"2"junk','CHR,BP\n1,2,3')) {
    input2 <- tempfile(tmpdir=d); writeBin(charToRaw(bad),input2); dest <- tempfile(tmpdir=d)
    expect_error(harmonise_file(input2,ref,dest,delimiter=',',batch_bytes=2),'quote|fields')
    expect_false(dir.exists(dest))
  }
  expect_error(harmonise_file(input,ref,file.path(d,'receipt-fail'),delimiter=',',acquisition_receipt=list(rows=3)), 'receipt')
  expect_false(dir.exists(file.path(d,'receipt-fail')))
})

test_that('typed statistical entities retain identity without small-variant effects', {
  d <- tempfile();dir.create(d);writeLines(c('>1','ACGT'),file.path(d,'f'))
  writeLines('chr\tpos\tref\talt\tid\taliases',file.path(d,'v'))
  ref <- prepare_reference(file.path(d,'f'),file.path(d,'v'),file.path(d,'r'),'GRCh37','toy')
  x <- data.frame(CHR='1',BP='2',A1='C',A2=c('<DEL>','C]2:5]','A','A','A'),
                  kind=c(NA,NA,'hla','repeat','gene_burden'),key=c('sv1','bnd1','HLA-A*02:01','HTT:CAG:40','GENE:mask:pLoF'),BETA='3')
  cols <- list(chr='CHR',pos='BP',other='A1',effect_allele='A2',effect='BETA',entity_type='kind',entity_id='key')
  y <- harmonise_sumstats(x,ref,columns=cols)
  expect_identical(y$source,x)
  expect_identical(y$normalized$entity_type,c('symbolic_sv','breakend','hla','repeat','gene_burden'))
  expect_identical(y$normalized$entity_id,x$key)
  expect_true(all(y$normalized$status=='specialist_preserved'))
  expect_true(all(is.na(y$normalized$effect)))
  expect_true(all(is.na(y$normalized$pos)))
})

test_that('specialist rows bypass small-variant candidate resolution', {
  d <- tempfile();dir.create(d);writeLines(c('>1','ACGT'),file.path(d,'f'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t2\tC\tA,G\trs2\t.'),file.path(d,'v'))
  ref <- prepare_reference(file.path(d,'f'),file.path(d,'v'),file.path(d,'r'),'GRCh37','toy')
  x <- data.frame(CHR='1',BP=2,SNP='rs2',A1='C',A2='A',kind='gene_burden',BETA=1)
  y <- harmonise_sumstats(x,ref,columns=list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',entity_type='kind'),candidate_limit=1)
  expect_identical(y$normalized$status,'specialist_preserved')
  expect_true(is.na(y$normalized$effect))
  expect_identical(y$normalized$entity_id,'rs2')
})

test_that('quoted decoder preserves UTF-8 and all field boundaries at arbitrary byte splits', {
  open_parser <- get('mss_parser_open_cpp',envir=environment(harmonise_file))
  feed_parser <- get('mss_parser_feed_cpp',envir=environment(harmonise_file))
  expected <- list(c('CHR','note','empty'),c('22','alpha,\tβ\r\nsaid "yes"',''),c('','"','tail'))
  for (delimiter in c('\t',',')) {
    encode <- function(row) paste(paste0('"',gsub('"','""',row,fixed=TRUE),'"'),collapse=delimiter)
    bytes <- charToRaw(paste(vapply(expected,encode,character(1)),collapse='\r\n'))
    for (size in c(1,2,7,31,65536)) {
      parser <- open_parser(delimiter,10000L);got <- list()
      for (i in seq.int(1,length(bytes),by=size)) got <- c(got,feed_parser(parser,bytes[i:min(length(bytes),i+size-1L)],FALSE))
      got <- c(got,feed_parser(parser,raw(),TRUE))
      expect_equal(got,expected)
    }
  }
  parser <- open_parser(',',8L)
  expect_error(feed_parser(parser,charToRaw('"123456789"'),TRUE),'max_record_bytes')
})
