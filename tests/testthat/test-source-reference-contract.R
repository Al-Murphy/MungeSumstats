test_that('declared source REF ALT strand prevents opposite indel event mapping', {
  d <- tempfile('source-contract-');dir.create(d)
  on.exit(unlink(d,recursive=TRUE),add=TRUE)
  writeLines(c('>1','CAAAAACGTACGT','>2','ACGT'),file.path(d,'reference.fa'))
  # Catalog contains an insertion; the source explicitly reports a deletion.
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t4\tA\tAA\trsIns\t.',
               '1\t10\tA\tC\trsSNP\t.','1\t8\tG\tA,T\trsMulti\t.',
               '2\t2\tC\tG\trsPal\t.'),file.path(d,'sites.tsv'))
  ref <- prepare_reference(file.path(d,'reference.fa'),file.path(d,'sites.tsv'),file.path(d,'index'),build='SYNTHETIC',version='trusted-source')
  columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',
                  source_ref='REF',source_alt='ALT',source_strand='STRAND')
  x <- data.frame(CHR='1',BP=4,SNP='unknown',A1='AA',A2='A',BETA=.125,FRQ=.25,REF='AA',ALT='A',STRAND='+')
  generic <- harmonise_sumstats(x,ref,frequency_type='EAF')$normalized
  expect_identical(generic$status,'ambiguous')
  expect_true(is.na(generic$effect))
  declared <- harmonise_sumstats(x,ref,columns=columns,frequency_type='EAF')
  trusted <- declared$normalized
  expect_identical(declared$provenance$implementation_version,'mss-preservation-4')
  expect_identical(declared$provenance$core_version,'mss-preservation-2')
  expect_identical(trusted$status,'sequence_only')
  expect_equal(trusted$pos,1)
  expect_identical(trusted$ref,'CA')
  expect_identical(trusted$alt,'C')
  expect_equal(trusted$effect,.125)
  expect_equal(trusted$eaf,.25)
  expect_true(is.na(trusted$id))
  x$A1 <- 'A';x$A2 <- 'AA';x$BETA <- -.125;x$FRQ <- .75
  unknown_site <- harmonise_sumstats(x,ref,columns=columns,frequency_type='EAF')$normalized
  expect_identical(unknown_site$status,'sequence_only')
  expect_equal(unknown_site$effect,.125)
  expect_true(is.na(unknown_site$eaf))
  x$A1 <- 'AA';x$A2 <- 'A';x$BETA <- .125;x$FRQ <- .25
  x$SNP <- 'rsIns'
  expect_identical(harmonise_sumstats(x,ref,columns=columns)$normalized$status,'source_id_conflict')
  # A known catalog identity still supports explicitly unordered allele input.
  known <- harmonise_sumstats(x,ref,frequency_type='EAF')$normalized
  expect_identical(known$status,'swap')
  expect_equal(known$effect,-.125)
  x <- data.frame(CHR='2',BP=2,SNP='rsPal',A1='C',A2='G',BETA=.125,FRQ=.25,REF='C',ALT='G',STRAND='+')
  pal <- harmonise_sumstats(x,ref,columns=columns,frequency_type='EAF')$normalized
  expect_identical(pal$status,'match')
  expect_equal(pal$effect,.125)
  x$STRAND <- '-';x$REF <- x$A1 <- 'G';x$ALT <- x$A2 <- 'C'
  reverse <- harmonise_sumstats(x,ref,columns=columns,frequency_type='EAF')$normalized
  expect_identical(reverse$status,'complement')
  expect_identical(reverse$ref,'C')
  multiallelic <- data.frame(CHR='1',BP=8,SNP='rsMulti',A1='T',A2='G',BETA=-.125,FRQ=.75,REF='G',ALT='T',STRAND='+')
  multi <- harmonise_sumstats(multiallelic,ref,columns=columns,frequency_type='EAF')$normalized
  expect_identical(multi$status,'swap')
  expect_equal(multi$effect,.125)
  expect_true(is.na(multi$eaf))
  x$STRAND <- 'unknown'
  expect_identical(harmonise_sumstats(x,ref,columns=columns)$normalized$status,'source_strand_unknown')
  x$STRAND <- '+';x$REF <- 'T';x$ALT <- 'A';x$A1 <- 'T';x$A2 <- 'A'
  mismatch <- harmonise_sumstats(x,ref,columns=columns)$normalized
  expect_identical(mismatch$status,'source_reference_mismatch')
  expect_true(is.na(mismatch$effect))
  expect_error(harmonise_sumstats(x,ref,columns=columns[names(columns)!='source_strand']),'together')
  x$REF <- 'C';x$ALT <- 'G';x$A1 <- 'A';x$A2 <- 'G'
  expect_identical(harmonise_sumstats(x,ref,columns=columns)$normalized$status,'source_allele_conflict')
  x$A1 <- 'C';x$A2 <- 'G';x$STRAND <- '+'
  write.table(x,file.path(d,'input.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  file <- harmonise_file(file.path(d,'input.tsv'),ref,file.path(d,'edition'),columns=columns,frequency_type='EAF')
  expect_equal(file$source_rows,1)
  expect_identical(read.delim(file.path(file$path,'normalized.tsv'))$status,'match')
  static <- columns[names(columns)!='source_strand']
  file_static <- harmonise_file(file.path(d,'input.tsv'),ref,file.path(d,'static-edition'),
                                columns=static,source_strand='+',frequency_type='EAF')
  expect_identical(file_static$source_strand,'+')
  expect_identical(file_static$implementation_version,'mss-preservation-4')
  expect_identical(readLines(file.path(d,'input.tsv')),readLines(file.path(file_static$path,'source.tsv')))
  metadata <- read.delim(file.path(file_static$path,'metadata.tsv'))
  expect_identical(metadata$value[metadata$key=='source_strand'],'+')
  expect_identical(metadata$value[metadata$key=='implementation_version'],'mss-preservation-4')
  expect_error(harmonise_sumstats(x,ref,columns=columns,source_strand='+'),'either')
})
