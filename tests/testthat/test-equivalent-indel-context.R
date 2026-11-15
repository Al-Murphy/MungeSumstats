test_that('equivalent indel anchors retain every third-allele context before EAF complement', {
  d <- tempfile('equivalent-context-');dir.create(d)
  on.exit(unlink(d,recursive=TRUE),add=TRUE)
  writeLines(c('>1','CAAAAAAACGT'),file.path(d,'reference.fa'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases',
               '1\t4\tAA\tA\trsDel\t.','1\t5\tAA\tA\trsDel\t.',
               '1\t5\tAA\tAG\trsThird\t.'),file.path(d,'sites.tsv'))
  ref <- prepare_reference(file.path(d,'reference.fa'),file.path(d,'sites.tsv'),file.path(d,'index'),
                           build='SYNTHETIC',version='equivalent-context')
  x <- data.frame(CHR='1',BP=4,SNP='rsDel',A1='A',A2='AA',BETA=-.125,FRQ=.75,REF='AA',ALT='A')
  generic <- harmonise_sumstats(x,ref,frequency_type='EAF')$normalized
  expect_identical(generic$status,'swap')
  expect_equal(generic$effect,.125)
  expect_true(is.na(generic$eaf))
  columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',source_ref='REF',source_alt='ALT')
  trusted <- harmonise_sumstats(x,ref,columns=columns,source_strand='+',frequency_type='EAF')$normalized
  expect_identical(trusted$status,'swap')
  expect_equal(trusted$effect,.125)
  expect_true(is.na(trusted$eaf))
  # Coordinate-only trusted normalization must also consult every equivalent
  # anchor, even when the original rsID is missing from the source.
  x$SNP <- NA_character_
  expect_true(is.na(harmonise_sumstats(x,ref,columns=columns,source_strand='+',frequency_type='EAF')$normalized$eaf))
})
