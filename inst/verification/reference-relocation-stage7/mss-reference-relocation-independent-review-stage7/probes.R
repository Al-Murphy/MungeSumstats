make_fixture <- function(format=2L, build='GRCh38') {
    root <- tempfile('independent-transfer-');dir.create(root)
    original <- file.path(root,'origin');dir.create(original)
    fasta <- file.path(original,'ref.fa')
    writeLines(c('>1','AACCGGTTAACCGGTT'),fasta,useBytes=TRUE)
    variants <- file.path(original,'variants.tsv')
    writeLines(c('chr\tpos\tref\talt\tid\taliases',
      '1\t3\tC\tA,G\trs30\told30;同義',
      '1\t9\tA\tAT\trs90\told90',
      '1\t13\tGG\tG\trs130\told130'),variants,useBytes=TRUE)
    reference <- prepare_reference(fasta,variants,file.path(original,'index'),
       build=build,version='independent-review',reference_format=paste0('mss-preservation-',format),
       contig_aliases=c(chr1='1',NC_000001.11='1'))
    input <- data.frame(CHR=c('chr1','1','NC_000001.11','chr1','1','chr1','chr1'),
      BP=c('0003','3','3','9','13',NA,'3'),
      SNP=c('rs30','old30','同義','old90','rs130','old30',NA),
      A1=c('C','A','C','A','GG','C','C'),A2=c('A','C','G','AT','G','A','A'),
      BETA=c('.2','-.2','.4','-.17','.03','.2','.2'),
      FRQ=c('.2','.8','.4','.1','.7','.2','.2'),
      raw='  untouched:000.00',stringsAsFactors=FALSE)
    list(root=root,original=original,fasta=normalizePath(fasta),reference=reference,input=input,
       hashes=hash_bundle(reference$path))
}
hash_bundle <- function(path) {
    files <- list.files(path,full.names=TRUE)
    stats::setNames(unname(tools::md5sum(files)),basename(files))
}
copy_fixture <- function(f) {
    copy <- file.path(f$root,'copy é 雪 with spaces');dir.create(copy)
    index <- file.path(copy,'index');dir.create(index)
    stopifnot(file.copy(f$fasta,file.path(copy,'different basename.fasta')),
      all(file.copy(list.files(f$reference$path,full.names=TRUE),index)))
    list(path=index,fasta=file.path(copy,'different basename.fasta'))
}
trim_manifest <- function(x) {x$path<-NULL;x$created_utc<-NULL;x}
expect_default_files <- function(a,b,projection=FALSE) {
    files <- c('normalized.tsv','source.tsv','source-artifact','audit.tsv','metadata.tsv',
      'columns.tsv','reference_checksums.tsv','reference_contigs.tsv','checksums.tsv',
      if(projection)c('projection.tsv','projection_audit.tsv'))
    expect_identical(unname(tools::md5sum(file.path(a$path,files))),
      unname(tools::md5sum(file.path(b$path,files))))
    expect_identical(trim_manifest(a),trim_manifest(b))
}

test_that('the original default API and all portable artifact fields remain exact', {
  expect_identical(names(formals(open_reference)),c('reference','fasta'))
  expect_null(formals(open_reference)$fasta)
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    handle<-open_reference(f$reference)
    expect_identical(open_reference(handle),handle)
    expect_identical(harmonise_sumstats(f$input,handle,frequency_type='EAF'),
      old$harmonise_sumstats(f$input,f$reference,frequency_type='EAF'))
    input<-file.path(f$root,'input.tsv');write.table(f$input,input,sep='\t',quote=TRUE,row.names=FALSE)
    a<-harmonise_file(input,f$reference,file.path(f$root,'new'),frequency_type='EAF',batch_size=2L)
    b<-old$harmonise_file(input,f$reference,file.path(f$root,'old'),frequency_type='EAF',batch_size=2L)
    expect_default_files(a,b)
    expect_identical(hash_bundle(f$reference$path),f$hashes)
  }
})

test_that('new physical streams preserve all selector and measured allele rows', {
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    before<-harmonise_sumstats(f$input,f$reference,frequency_type='EAF')
    first<-open_reference(f$reference)
    moved<-copy_fixture(f)
    # Normal copies and different basenames require no symlink privileges.
    second<-open_reference(f$reference,fasta=moved$fasta)
    expect_false(identical(second$.native,first$.native))
    expect_identical(second$path,f$reference$path)
    expect_identical(second$fasta,normalizePath(moved$fasta))
    expect_identical(second$checksums,first$checksums)
    expect_identical(open_reference(second),second)
    expect_identical(harmonise_sumstats(f$input,second,frequency_type='EAF')$normalized,before$normalized)
    expect_identical(before$normalized$effect[1:2],c(.2,.2))
    expect_identical(before$source,f$input)
    # Force new index path explicitly; no old-path inference is permitted.
    rm(first);invisible(gc());unlink(f$original,recursive=TRUE)
    expect_false(file.exists(f$fasta))
    expect_error(open_reference(moved$path),'checksum|size|exist')
    handle<-open_reference(moved$path,fasta=moved$fasta)
    after<-harmonise_sumstats(f$input,handle,frequency_type='EAF')
    expect_identical(after$normalized,before$normalized)
    expect_identical(after$source,f$input)
    expect_identical(after$provenance$reference_checksums,before$provenance$reference_checksums)
    expect_identical(hash_bundle(moved$path),f$hashes)
  }
})

test_that('explicit rebinding distrusts poisoned and serialized handles and stale metadata', {
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    handle<-open_reference(f$reference);moved<-copy_fixture(f)
    poisoned<-handle;poisoned$checksums[]<-'poison';poisoned$sizes[]<-0
    poisoned$files<-'nonexistent';poisoned$fasta<-'nonexistent';poisoned$build<-'wrong'
    poisoned$contig_aliases<-c(chr1='wrong');poisoned$fasta_binding<-list(runtime_path='wrong')
    poisoned$.file_signature<-NULL
    rebound<-open_reference(poisoned,fasta=moved$fasta)
    expect_identical(rebound$checksums,handle$checksums)
    expect_identical(rebound$build,handle$build)
    expect_identical(rebound$contig_aliases,handle$contig_aliases)
    expect_false(identical(rebound$.native,handle$.native))
    restored<-unserialize(serialize(rebound,NULL))
    # A dead external pointer cannot be reused directly.
    expect_error(harmonise_sumstats(f$input,restored),'no longer live|external pointer')
    reopened<-open_reference(restored,fasta=moved$fasta)
    expect_false(identical(reopened$.native,rebound$.native))
    expect_identical(harmonise_sumstats(f$input,reopened)$normalized,
      harmonise_sumstats(f$input,handle)$normalized)
    # Alter only the real manifest after the handle was opened.
    manifest<-file.path(f$reference$path,'manifest.rds');value<-readRDS(manifest)
    value$checksums[['fasta']]<-paste(rep('0',32),collapse='');saveRDS(value,manifest)
    expect_error(open_reference(handle,fasta=moved$fasta),'checksum')
    value$core_version<-'mss-preservation-999';saveRDS(value,manifest)
    expect_error(open_reference(handle,fasta=moved$fasta),'Unsupported reference schema')
  }
})

test_that('every binding refuses corrupt bytes including same-size alias and FASTA mutations', {
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    moved<-copy_fixture(f)
    for(name in c('fasta.fai','keys','offsets','records','source.tsv','contig_aliases.tsv','core-version')) {
      path<-file.path(moved$path,name);bytes<-readBin(path,'raw',file.info(path)$size)
      original<-bytes;bytes[[1]]<-as.raw(bitwXor(as.integer(bytes[[1]]),1L));writeBin(bytes,path)
      # Corrupting a quoted alias header also emits the existing read.delim
      # warning before it refuses; check refusal without treating it as a new warning.
      expect_error(suppressWarnings(open_reference(moved$path,fasta=moved$fasta)),'checksum|marker|alias')
      writeBin(original,path)
    }
    bytes<-readBin(moved$fasta,'raw',file.info(moved$fasta)$size);bytes[[4]]<-charToRaw('T')
    writeBin(bytes,moved$fasta)
    expect_equal(file.info(moved$fasta)$size,file.info(f$fasta)$size)
    expect_error(open_reference(moved$path,fasta=moved$fasta),'checksum')
    file.copy(f$fasta,moved$fasta,overwrite=TRUE)
    expect_s3_class(open_reference(moved$path,fasta=moved$fasta),'mss_reference_handle')
    expect_identical(hash_bundle(f$reference$path),f$hashes)
  }
})

test_that('malformed FAI syntax is refused and the inherited semantic-width limitation is explicit', {
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    moved<-copy_fixture(f);fai<-file.path(moved$path,'fasta.fai')
    writeLines('1\t16\t3\t0',fai)
    manifest<-file.path(moved$path,'manifest.rds');value<-readRDS(manifest)
    value$checksums[['fasta.fai']]<-unname(tools::md5sum(fai))
    value$sizes[['fasta.fai']]<-file.info(fai)$size;saveRDS(value,manifest)
    expect_error(open_reference(moved$path,fasta=moved$fasta),'FAI|FASTA|index|line')
    # A trusted manifest replacement is outside relocation's integrity contract.
    # Establish that semantic zero-width FAI acceptance predates this R change.
    writeLines('1\t16\t3\t0\t0',fai)
    value$checksums[['fasta.fai']]<-unname(tools::md5sum(fai))
    value$sizes[['fasta.fai']]<-file.info(fai)$size;saveRDS(value,manifest)
    expect_s3_class(open_reference(moved$path,fasta=moved$fasta),'mss_reference_handle')
    expect_s3_class(old$open_reference(moved$path),'mss_reference_handle')
  }
})

test_that('runtime changes and both verification-open interval edges invalidate explicit bindings', {
  native_open<-mss_open_reference_cpp;verify_original<-.mss_verify_reference
  on.exit(assign('mss_open_reference_cpp',native_open,envir=.GlobalEnv),add=TRUE)
  on.exit(assign('.mss_verify_reference',verify_original,envir=.GlobalEnv),add=TRUE)
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    moved<-copy_fixture(f);handle<-open_reference(moved$path,fasta=moved$fasta)
    Sys.setFileTime(moved$fasta,Sys.time()+60)
    expect_error(harmonise_sumstats(f$input,handle,verify=FALSE),'checksum/stat')
    file.copy(f$fasta,moved$fasta,overwrite=TRUE)
    assign('.mss_verify_reference',function(reference,verify) {
      answer<-verify_original(reference,verify)
      Sys.setFileTime(reference$fasta,Sys.time()+120)
      answer
    },envir=.GlobalEnv)
    expect_error(open_reference(moved$path,fasta=moved$fasta),'changed while opening')
    assign('.mss_verify_reference',verify_original,envir=.GlobalEnv)
    assign('mss_open_reference_cpp',function(path,fasta) {
      answer<-native_open(path,fasta);Sys.setFileTime(fasta,Sys.time()+180);answer
    },envir=.GlobalEnv)
    expect_error(open_reference(moved$path,fasta=moved$fasta),'changed while opening')
    assign('mss_open_reference_cpp',native_open,envir=.GlobalEnv)
  }
})

test_that('source-only and target-only projection bindings retain old scientific receipt fields', {
  for(format in 2:3) {
    s<-make_fixture(format,'GRCh37');t<-make_fixture(5L-format,'GRCh38')
    on.exit(unlink(s$root,recursive=TRUE),add=TRUE);on.exit(unlink(t$root,recursive=TRUE),add=TRUE)
    chain<-file.path(s$root,'identity.chain')
    writeLines(c('chain 100 1 16 + 0 16 1 16 + 0 16 1','16'),chain)
    baseline<-prepare_projection(chain,s$reference,t$reference,'GRCh37','GRCh38')
    baseline_old<-old$prepare_projection(chain,s$reference,t$reference,'GRCh37','GRCh38')
    expect_identical(baseline$receipt,baseline_old$receipt)
    before<-harmonise_sumstats(s$input,s$reference,projection=baseline,frequency_type='EAF')
    input<-file.path(s$root,'input.tsv');write.table(s$input,input,sep='\t',quote=TRUE,row.names=FALSE)
    edition<-harmonise_file(input,s$reference,file.path(s$root,'new-projection'),projection=baseline,frequency_type='EAF')
    old_edition<-old$harmonise_file(input,s$reference,file.path(s$root,'old-projection'),projection=baseline_old,frequency_type='EAF')
    expect_default_files(edition,old_edition,TRUE)
    ms<-copy_fixture(s);mt<-copy_fixture(t)
    hs<-open_reference(ms$path,fasta=ms$fasta);ht<-open_reference(mt$path,fasta=mt$fasta)
    for(side in c('source','target')) {
      sr<-if(side=='source')hs else s$reference
      tr<-if(side=='target')ht else t$reference
      p<-prepare_projection(chain,sr,tr,'GRCh37','GRCh38')
      expect_identical(p$receipt[names(baseline$receipt)],baseline$receipt)
      expect_identical(is.null(p$receipt$source_fasta_binding),side!='source')
      expect_identical(is.null(p$receipt$target_fasta_binding),side!='target')
      after<-harmonise_sumstats(s$input,sr,projection=p,frequency_type='EAF')
      expect_identical(after$source,before$source)
      expect_identical(after$normalized,before$normalized)
      expect_identical(after$provenance$reference_checksums,before$provenance$reference_checksums)
      # A changed projection side must be rejected even if harmonisation uses the other original side.
      chosen<-if(side=='source')ms$fasta else mt$fasta
      when<-file.info(chosen)$mtime;Sys.setFileTime(chosen,Sys.time()+240)
      expect_error(harmonise_sumstats(s$input,sr,projection=p),'checksum/stat')
      Sys.setFileTime(chosen,when)
    }
    expect_identical(hash_bundle(ms$path),s$hashes)
    expect_identical(hash_bundle(mt$path),t$hashes)
  }
})

test_that('empty-row streams retain explicit operational receipt and immutable scientific identity', {
  for(format in 2:3) {
    f<-make_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    moved<-copy_fixture(f);handle<-open_reference(moved$path,fasta=moved$fasta)
    empty<-f$input[FALSE,,drop=FALSE]
    input<-file.path(f$root,'empty.tsv');write.table(empty,input,sep='\t',quote=TRUE,row.names=FALSE)
    result<-harmonise_file(input,handle,file.path(f$root,'empty-out'))
    expect_identical(result$source_rows,0)
    expect_identical(result$reference_fasta_binding,handle$fasta_binding)
    expect_identical(result$reference_checksums,f$reference$checksums)
    expect_equal(nrow(read.delim(file.path(result$path,'normalized.tsv'))),0)
    expect_identical(hash_bundle(moved$path),f$hashes)
  }
})
