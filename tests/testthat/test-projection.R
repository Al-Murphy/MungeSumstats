projection_fixture <- function() {
  d <- tempfile();dir.create(d)
  make <- function(name, seq, sites) {
    writeLines(c('>1',seq),file.path(d,paste0(name,'.fa')))
    writeLines(c('chr\tpos\tref\talt\tid\taliases',sites),file.path(d,paste0(name,'.tsv')))
    prepare_reference(file.path(d,paste0(name,'.fa')),file.path(d,paste0(name,'.tsv')),file.path(d,name),if(name=='src') 'GRCh37' else 'GRCh38','toy')
  }
  src <- make('src','CAAAAACGTACGTACGTACG',c('1\t7\tC\tA,G\trs7\t.','1\t4\tAA\tA\trsDel\t.'))
  target <- make('target','CAAAAACGTATACGCGTACG',c('1\t7\tC\tA,G\trs7\t.','1\t4\tAA\tA\trsDel\t.'))
  chain <- file.path(d,'toy.chain')
  writeLines(c('chain 100 1 20 + 0 10 1 20 + 0 10 1','10','',
               'chain 100 1 20 + 10 14 1 20 - 6 10 2','4','',
               'chain 100 1 20 + 14 17 1 20 + 14 17 3','1 1 1','1','',
               'chain 100 1 20 + 18 19 1 20 + 18 19 4','1','',
               'chain 100 1 20 + 18 19 1 20 + 10 11 5','1'),chain)
  list(d=d,source=src,target=target,chain=chain)
}
test_that('projection preserves native identity, signed effects and nullable outcomes', {
  f <- projection_fixture()
  expect_true(exists('prepare_projection',mode='function'))
  if(!exists('prepare_projection',mode='function')) return(invisible(NULL))
  p <- prepare_projection(f$chain,f$source,f$target,'GRCh37','GRCh38')
  x <- data.frame(CHR='1',BP=c('7','4','11','15','18','19'),REF=c('C','AA','C','CG','A','C'),ALT=c('A','A','A','C','C','A'),A1=c('A','AA','C','CG','A','C'),A2=c('C','A','A','C','C','A'),BETA='-.125',FRQ='.75')
  cols <- list(chr='CHR',pos='BP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',source_ref='REF',source_alt='ALT')
  a <- harmonise_sumstats(x,f$source,columns=cols,source_strand='+',frequency_type='EAF')
  b <- harmonise_sumstats(x,f$source,columns=cols,source_strand='+',frequency_type='EAF',projection=p)
  expect_identical(b$source,x)
  expect_identical(b$normalized[names(a$normalized)[!grepl('grch38',names(a$normalized))]],a$normalized[names(a$normalized)[!grepl('grch38',names(a$normalized))]])
  expect_identical(b$normalized$projection_status,c('projected','projected','projected','chain_gap','unmapped','multimap'))
  expect_equal(b$normalized$projected_pos[1:3],c(7,1,14))
  expect_identical(b$normalized$projected_ref[1:3],c('C','CA','G'))
  expect_identical(b$normalized$projected_alt[1:3],c('A','C','T'))
  expect_equal(b$normalized$projected_effect[1:3],c(.125,-.125,-.125))
  expect_true(is.na(b$normalized$projected_eaf[1]))
  expect_equal(b$normalized$projected_eaf[2:3],c(.75,.75))
  expect_true(all(is.na(b$normalized$projected_pos[4:6])))
  expect_equal(b$normalized$grch38_pos,b$normalized$projected_pos)
  expect_error(prepare_projection(f$chain,f$source,f$target,'wrong','GRCh38'),'build')
  cat('\n',file=f$chain,append=TRUE)
  expect_error(harmonise_sumstats(x,f$source,columns=cols,source_strand='+',projection=p),'changed|checksum')
})

test_that('negative indel spans left align and target third alleles suppress complemented EAF', {
  f <- projection_fixture()
  # Whole contig reverse complement: deletion in the five-A run becomes a
  # deletion in the five-T run at target 14, GT>G, after left alignment.
  writeLines(c('>1','CGTACGTACGTACGTTTTTG'),file.path(f$d,'reverse.fa'))
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t14\tGT\tG\trsDel\t.'),file.path(f$d,'reverse.tsv'))
  target <- prepare_reference(file.path(f$d,'reverse.fa'),file.path(f$d,'reverse.tsv'),file.path(f$d,'reverse'),'GRCh38','reverse')
  chain <- file.path(f$d,'reverse.chain')
  writeLines(c('chain 100 1 20 + 0 20 1 20 - 0 20 1','20'),chain)
  p <- prepare_projection(chain,f$source,target,'GRCh37','GRCh38')
  x <- data.frame(CHR='1',BP=4,REF='AA',ALT='A',A1='AA',A2='A',BETA=.125,FRQ=.2)
  cols <- list(chr='CHR',pos='BP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',source_ref='REF',source_alt='ALT')
  y <- harmonise_sumstats(x,f$source,columns=cols,source_strand='+',frequency_type='EAF',projection=p)$normalized
  expect_equal(y$projected_pos,14)
  expect_identical(y$projected_ref,'GT');expect_identical(y$projected_alt,'G')
  expect_equal(y$projected_effect,.125);expect_equal(y$projected_eaf,.2)
  # Source says biallelic but target has a third allele at the equivalent site.
  writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t7\tC\tA\trs7\t.'),file.path(f$d,'biallelic.tsv'))
  source <- prepare_reference(f$source$fasta,file.path(f$d,'biallelic.tsv'),file.path(f$d,'biallelic'),'GRCh37','biallelic')
  p <- prepare_projection(f$chain,source,f$target,'GRCh37','GRCh38')
  x <- data.frame(CHR='1',BP=7,REF='C',ALT='A',A1='A',A2='C',BETA=.125,FRQ=.2)
  y <- harmonise_sumstats(x,source,columns=cols,source_strand='+',frequency_type='EAF',projection=p)$normalized
  expect_equal(y$effect,-.125);expect_equal(y$eaf,.8)
  expect_equal(y$projected_effect,-.125);expect_true(is.na(y$projected_eaf))
})

test_that('chain gaps, reference conflicts and malformed build assets stay explicit', {
  f <- projection_fixture()
  chain <- file.path(f$d,'gap.chain')
  writeLines(c('chain 100 1 20 + 0 20 1 20 + 0 20 1','10 4 4','6'),chain)
  p <- prepare_projection(chain,f$source,f$target,'GRCh37','GRCh38')
  x <- data.frame(CHR='1',BP=c(11,14,1),REF=c('C','A','A'),ALT=c('A','G','C'),A1=c('C','A','A'),A2=c('A','G','C'))
  cols <- list(chr='CHR',pos='BP',other='A1',effect_allele='A2',source_ref='REF',source_alt='ALT')
  y <- harmonise_sumstats(x,f$source,columns=cols,source_strand='+',projection=p)$normalized
  expect_identical(y$projection_status,c('chain_gap','chain_gap','native_unresolved'))
  writeLines(c('chain 100 1 20 + 0 20 1 20 + 0 20 1','20'),chain)
  p <- prepare_projection(chain,f$source,f$target,'GRCh37','GRCh38')
  y <- harmonise_sumstats(x,f$source,columns=cols,source_strand='+',projection=p)$normalized
  expect_identical(y$projection_status,c('target_ref_conflict','target_ref_conflict','native_unresolved'))
  expect_true(all(is.na(y$projected_effect)))
  writeLines(c('chain 100 1 21 + 0 20 1 20 + 0 20 1','20'),chain)
  expect_error(prepare_projection(chain,f$source,f$target,'GRCh37','GRCh38'),'length|size')
  writeLines(c('chain 100 1 20 + 0 20 1 20 + 0 20 1','10 4 4'),chain)
  expect_error(prepare_projection(chain,f$source,f$target,'GRCh37','GRCh38'),'Truncated')
})

test_that('chain IDs may repeat on distinct source contigs in bundled Ensembl chains', {
  f <- projection_fixture()
  writeLines(c('chain 100 1 20 + 0 20 1 20 + 0 20 1','20','',
               'chain 100 2 20 + 0 20 1 20 + 0 20 1','20'),f$chain)
  expect_s3_class(prepare_projection(f$chain,f$source,f$target,'GRCh37','GRCh38'),'mss_projection')
})

test_that('explicit aliases cannot collapse distinct chains into a false unique indel mapping', {
  f <- projection_fixture()
  writeLines(c('chain 100 1 20 + 0 20 1 20 + 0 20 1','20','',
               'chain 100 chr1 20 + 0 20 1 20 + 0 20 1','1 1 1','18'),f$chain)
  p <- prepare_projection(f$chain,f$source,f$target,'GRCh37','GRCh38',source_contigs=c(chr1='1'))
  x <- data.frame(CHR='1',BP=4,REF='AA',ALT='A',A1='AA',A2='A')
  cols <- list(chr='CHR',pos='BP',other='A1',effect_allele='A2',source_ref='REF',source_alt='ALT')
  y <- harmonise_sumstats(x,f$source,columns=cols,source_strand='+',projection=p)$normalized
  expect_identical(y$projection_status,'multimap')
  expect_true(is.na(y$projected_pos))
})
