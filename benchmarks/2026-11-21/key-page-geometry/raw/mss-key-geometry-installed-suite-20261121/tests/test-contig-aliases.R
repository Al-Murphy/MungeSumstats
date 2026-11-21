test_that('explicit contig aliases preserve source labels and resolve reference labels', {
  d <- tempfile('contig-alias-');dir.create(d)
  on.exit(unlink(d,recursive=TRUE),add=TRUE)
  writeLines(c('>1','ACGTACGT','>2','ACGTACGT'),file.path(d,'ref.fa'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t2\tC\tA\trs2\t.'),file.path(d,'sites.tsv'))
  aliases <- c('1'='1',chr1='1',NC_000001.10='1')
  ref <- prepare_reference(file.path(d,'ref.fa'),file.path(d,'sites.tsv'),
                           file.path(d,'index'),build='GRCh37',version='explicit-report',
                           contig_aliases=aliases)
  x <- data.frame(CHR=c('1','chr1','NC_000001.10','unknown'),BP=2,SNP=NA_character_,
                  A1='C',A2='A',BETA=.25,FRQ=.2)
  result <- harmonise_sumstats(x,ref,frequency_type='EAF')
  expect_identical(result$source,x)
  expect_identical(result$normalized$chr[1:3],rep('1',3))
  expect_identical(result$normalized$id[1:3],rep('rs2',3))
  expect_equal(result$normalized$effect[1:3],rep(.25,3))
  expect_identical(result$normalized$status[4],'unmapped')
  expect_identical(result$provenance$contig_mapping_implementation,'mss-contigs-1')
  expect_true('contig_aliases.tsv' %in% names(ref$checksums))
  handle <- open_reference(ref)
  expect_identical(harmonise_sumstats(x,handle)$normalized$chr,result$normalized$chr)
  input <- file.path(d,'source.tsv');write.table(x,input,sep='\t',row.names=FALSE,quote=FALSE)
  edition <- harmonise_file(input,handle,file.path(d,'edition'),frequency_type='EAF')
  expect_identical(readBin(input,'raw',n=file.info(input)$size),
                   readBin(file.path(d,'edition','source-artifact'),'raw',n=file.info(input)$size))
  expect_true('reference_contigs.tsv' %in% names(edition$checksums))
  table <- read.delim(file.path(d,'edition','reference_contigs.tsv'),colClasses='character')
  expect_setequal(table$alias,names(aliases))
  expect_identical(edition$contig_mapping_implementation,'mss-contigs-1')
  writeLines('changed',file.path(ref$path,'contig_aliases.tsv'))
  expect_error(harmonise_sumstats(x,handle),'checksum/stat')
})

test_that('contig mappings cannot shadow canonical names or invent destinations', {
  d <- tempfile('invalid-contig-alias-');dir.create(d)
  on.exit(unlink(d,recursive=TRUE),add=TRUE)
  writeLines(c('>1','ACGT','>2','ACGT'),file.path(d,'ref.fa'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t2\tC\tA\trs2\t.'),file.path(d,'sites.tsv'))
  for (aliases in list(c('1'='2'),c(chr1='missing'),setNames(c('1','2'),c('same','same')),
                       c('1'),c(chr1=NA_character_))) {
    path <- tempfile('index-',tmpdir=d)
    expect_error(prepare_reference(file.path(d,'ref.fa'),file.path(d,'sites.tsv'),path,
                                   build='B',version='V',contig_aliases=aliases),'[Cc]ontig|[Aa]lias')
    expect_false(dir.exists(path))
  }
  old <- prepare_reference(file.path(d,'ref.fa'),file.path(d,'sites.tsv'),file.path(d,'old'),build='B',version='V')
  x <- data.frame(CHR='chr1',BP=2,SNP=NA_character_,A1='C',A2='A',BETA=.25)
  expect_identical(harmonise_sumstats(x,old)$normalized$status,'unmapped')
  expect_false('contig_aliases.tsv' %in% names(old$checksums))
})
