fixture <- function(format) {
    root <- tempfile('independent-fai-');dir.create(root)
    fasta <- file.path(root,'reference.fa')
    writeBin(charToRaw('>empty\n>1\nCAAAAACGTACGT\n'),fasta)
    sites <- file.path(root,'sites.tsv')
    writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t8\tG\tA,T\trs8\told8',
                 '1\t4\tAA\tA\tdel\told_del'),sites)
    ref <- prepare_reference(fasta,sites,file.path(root,'index'),build='GRCh38',version='independent-fai',
      reference_format=paste0('mss-preservation-',format),contig_aliases=c(chr1='1'))
    x <- data.frame(CHR='chr1',BP=c('8','8','4','8'),SNP=c('old8','rs8','old_del','rs8'),
      A1=c('G','A','AA','A'),A2=c('A','G','A','T'),BETA=c('.3','-.3','.125','.8'),
      FRQ=c('.2','.8','.25','.7'),raw='raw:000.00')
    list(root=root,fasta=fasta,ref=ref,input=x)
}
repin <- function(f,bytes,name) {
    path<-file.path(f$root,name);dir.create(path)
    stopifnot(all(file.copy(list.files(f$ref$path,full.names=TRUE),path)))
    file<-file.path(path,'fasta.fai');writeBin(bytes,file)
    manifest_path<-file.path(path,'manifest.rds');m<-readRDS(manifest_path)
    m$checksums[['fasta.fai']]<-unname(tools::md5sum(file));m$sizes[['fasta.fai']]<-file.info(file)$size
    saveRDS(m,manifest_path,version=2);path
}
bundle_hashes<-function(path) {
    paths<-list.files(path,full.names=TRUE)
    stats::setNames(unname(tools::md5sum(paths)),basename(paths))
}

test_that('valid retained handles preserve complete baseline results and scientific editions', {
  for(format in 2:3) {
    f<-fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    hashes<-bundle_hashes(f$ref$path)
    candidate<-open_reference(f$ref);baseline<-old$open_reference(f$ref)
    for(batch in c(1L,4L)) {
      actual<-harmonise_sumstats(f$input,candidate,frequency_type='EAF',batch_size=batch)
      expected<-old$harmonise_sumstats(f$input,baseline,frequency_type='EAF',batch_size=batch)
      expect_identical(actual,expected)
      expect_identical(actual$source,f$input)
      expect_identical(ncol(actual$normalized),19L)
      expect_identical(actual$normalized$effect[1:3],c(.3,.3,.125))
      expect_identical(actual$normalized$status[[4]],'alt_alt_contrast')
    }
    input<-file.path(f$root,'input.tsv');write.table(f$input,input,sep='\t',quote=TRUE,row.names=FALSE)
    after<-harmonise_file(input,candidate,file.path(f$root,'current'),frequency_type='EAF',batch_size=1L)
    before<-old$harmonise_file(input,baseline,file.path(f$root,'baseline'),frequency_type='EAF',batch_size=1L)
    files<-c('normalized.tsv','source.tsv','source-artifact','audit.tsv','metadata.tsv','columns.tsv',
             'reference_checksums.tsv','reference_contigs.tsv','checksums.tsv')
    expect_identical(unname(tools::md5sum(file.path(after$path,files))),
                     unname(tools::md5sum(file.path(before$path,files))))
    after$created_utc<-NULL;after$path<-NULL;before$created_utc<-NULL;before$path<-NULL
    expect_identical(after,before)
    expect_identical(bundle_hashes(f$ref$path),hashes)
    expect_identical(open_reference(candidate),candidate)
    rm(candidate,baseline);invisible(gc())
  }
})

test_that('fresh and rebound public handles refuse repinned unsafe geometry after prior caching', {
  for(format in 2:3) {
    f<-fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    prior<-open_reference(f$ref)
    actual<-harmonise_sumstats(f$input,prior)
    body<-c('1\t1\t3\t1\t18446744073709551616\n',
      '1\t18446744073709551615\t3\t18446744073709551615\t18446744073709551615\n',
      '1\t2\t3\t1\t18446744073709551615\n','1\t1\t18446744073709551614\t1\t1\n',
      '1\t0\t0\t0\t18446744073709551615\n','1\t1\t3\t+1\t2\n',
      '1\t1\t3\t1 \t2\n','1\t1\t3\t1\t2\r\r\n')
    for(i in seq_along(body)) {
      path<-repin(f,charToRaw(body[[i]]),paste0('unsafe-',i))
      expect_error(open_reference(path),'FASTA index')
      expect_error(open_reference(path,fasta=f$fasta),'FASTA index')
      poisoned<-prior;poisoned$path<-path
      expect_error(open_reference(poisoned,fasta=f$fasta),'FASTA index')
    }
    # A previously verified handle cannot hide an on-disk FAI mutation.
    old_manifest<-file.path(f$ref$path,'manifest.rds');m<-readRDS(old_manifest)
    index<-file.path(f$ref$path,'fasta.fai');writeLines('1\t13\t10\t0\t0',index)
    m$checksums[['fasta.fai']]<-unname(tools::md5sum(index));m$sizes[['fasta.fai']]<-file.info(index)$size
    saveRDS(m,old_manifest,version=2)
    expect_error(harmonise_sumstats(f$input,prior,verify=FALSE),'checksum/stat')
    expect_error(open_reference(prior,fasta=f$fasta),'FASTA index')
    rm(prior);invisible(gc())
  }
})

test_that('maximum unused stride and empty EOF offsets remain safe supported geometry', {
  for(format in 2:3) {
    f<-fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    # Only one addressed row: the declared stride is never used to reach a second line.
    length<-file.info(f$fasta)$size
    bytes<-charToRaw(paste0('1\t13\t10\t13\t18446744073709551615\nempty\t0\t',length,'\t0\t0'))
    path<-repin(f,bytes,'max-unused-stride')
    handle<-open_reference(path,fasta=f$fasta)
    baseline<-open_reference(f$ref)
    expect_identical(harmonise_sumstats(f$input,handle)$normalized,harmonise_sumstats(f$input,baseline)$normalized)
    expect_identical(handle$checksums[['fasta.fai']],unname(tools::md5sum(file.path(path,'fasta.fai'))))
    expect_identical(handle$fasta_binding$verification,'full_reference_md5')
    rm(handle,baseline);invisible(gc())
  }
})
