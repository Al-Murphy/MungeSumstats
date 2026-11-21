.fai_safety_fixture <- function(format, bytes=charToRaw('>one\nACGT\nAC\n>two\nT\n')) {
    root <- tempfile('fai-numeric-');dir.create(root)
    fasta <- file.path(root,'reference.fa');writeBin(bytes,fasta)
    sites <- file.path(root,'sites.tsv');writeBin(charToRaw('chr\tpos\tref\talt\tid\taliases\n'),sites)
    ref <- prepare_reference(fasta,sites,file.path(root,'reference'),build='SYNTHETIC',version='fai-safety',
                             reference_format=paste0('mss-preservation-',format))
    list(root=root,fasta=fasta,reference=ref)
}
.fai_safety_repin <- function(fixture,bytes,name) {
    path <- file.path(fixture$root,name);dir.create(path)
    stopifnot(all(file.copy(list.files(fixture$reference$path,full.names=TRUE),path)))
    file <- file.path(path,'fasta.fai');writeBin(bytes,file)
    manifest <- readRDS(file.path(path,'manifest.rds'))
    manifest$checksums[['fasta.fai']] <- unname(tools::md5sum(file))
    manifest$sizes[['fasta.fai']] <- unname(file.info(file)$size)
    saveRDS(manifest,file.path(path,'manifest.rds'),version=2)
    path
}

test_that('re-pinned malformed FAI rows are refused by the native gate in both formats', {
    valid <- 'one\t6\t5\t4\t5\ntwo\t1\t18\t1\t2\n'
    bad <- list(zero_bases='one\t6\t5\t0\t5\n',zero_bytes='one\t6\t5\t4\t0\n',
                narrow_bytes='one\t6\t5\t4\t3\n',bases_exceed_size='one\t6\t5\t7\t8\n',
                offset_beyond_file='one\t6\t999\t4\t5\n',last_base_beyond_file='one\t100\t5\t4\t5\n',
                offset_overflow='one\t6\t18446744073709551615\t4\t5\n',
                multiplication_overflow='one\t18446744073709551615\t5\t1\t18446744073709551615\n',
                duplicate_name=paste0(valid,'one\t6\t5\t4\t5\n'),empty_name='\t6\t5\t4\t5\n',
                whitespace_name='one bad\t6\t5\t4\t5\n',empty_index='',empty_numeric='one\t\t5\t4\t5\n',
                negative='one\t-1\t5\t4\t5\n',plus='one\t+6\t5\t4\t5\n',suffix='one\t6x\t5\t4\t5\n',
                leading_space='one\t 6\t5\t4\t5\n',uint_overflow='one\t18446744073709551616\t5\t4\t5\n',
                extra_field='one\t6\t5\t4\t5\t0\n',inconsistent_empty='empty\t0\t0\t1\t1\n',
                trailing_space='one\t6\t5\t4\t5 \n',hex_suffix='one\t6\t0x5\t4\t5\n',
                empty_offset_beyond_file='empty\t0\t999\t0\t0\n',
                NUL_name=c(charToRaw('on'),as.raw(0),charToRaw('e\t6\t5\t4\t5\n')))
    for(format in 2:3) {
        f <- .fai_safety_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        for(name in names(bad)) {
            bytes <- bad[[name]];if(is.character(bytes))bytes <- charToRaw(bytes)
            path <- .fai_safety_repin(f,bytes,name)
            manifest <- readRDS(file.path(path,'manifest.rds'))
            expect_identical(unname(tools::md5sum(file.path(path,'fasta.fai'))),manifest$checksums[['fasta.fai']],info=name)
            handle <- tryCatch(open_reference(path),error=function(e)e)
            expect_true(inherits(handle,'error'),info=paste(format,name,'must refuse before sequence arithmetic'))
            if(inherits(handle,'error'))
                expect_match(conditionMessage(handle),'FASTA index',info=paste(format,name,'native refusal, not checksum gate'))
            rm(handle);invisible(gc())
        }
    }
})

test_that('valid indexed LF CRLF and unterminated FASTAs retain all public scientific fields', {
    controls <- list(LF=charToRaw('>one\nACGT\nAC\n'),CRLF=charToRaw('>one\r\nACGT\r\nAC\r\n'),
                     no_final_newline=charToRaw('>one\nACGT\nAC'),single_unterminated=charToRaw('>one\nACGTAC'),
                     empty_contig=charToRaw('>empty\n>one\nACGT\nAC\n'),UTF8=charToRaw('>雪\nACGT\nAC\n'))
    for(format in 2:3)for(name in names(controls)) {
        f <- .fai_safety_fixture(format,controls[[name]]);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        chr <- if(name=='UTF8') '雪' else 'one'
        ref <- c('A','C','G','T','A','C');alt <- c('C','A','A','C','C','A')
        x <- data.frame(CHR=chr,BP=1:6,A1=ref,A2=alt,BETA=seq_len(6)/100,FRQ=.2,untouched='raw:0001')
        columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',
                        source_ref='A1',source_alt='A2')
        handle <- open_reference(f$reference)
        result <- harmonise_sumstats(x,handle,columns=columns,source_strand='+',frequency_type='EAF')
        expect_identical(result$source,x)
        expect_identical(ncol(result$normalized),19L)
        expect_identical(result$normalized$status,rep('sequence_only',6))
        expect_equal(result$normalized$pos,1:6)
        expect_identical(result$normalized$ref,ref);expect_identical(result$normalized$alt,alt)
        expect_equal(result$normalized$effect,x$BETA);expect_equal(result$normalized$eaf,x$FRQ)
        expect_identical(result$provenance$implementation_version,'mss-preservation-7')
        # Leading zeros and CRLF FAI text are accepted decimal representations.
        fai <- readLines(file.path(f$reference$path,'fasta.fai'),warn=FALSE)
        padded <- vapply(strsplit(fai,'\t',fixed=TRUE),function(v)paste(c(v[[1]],paste0('000',v[-1])),collapse='\t'),character(1))
        path <- .fai_safety_repin(f,charToRaw(paste0(paste(padded,collapse='\r\n'),'\r\n')),'leading-zero-index')
        leading <- open_reference(path,fasta=f$fasta)
        output <- harmonise_sumstats(x,leading,columns=columns,source_strand='+',frequency_type='EAF')
        expect_identical(output$source,result$source)
        expect_identical(output$normalized,result$normalized)
        expect_identical(output$provenance$implementation_version,result$provenance$implementation_version)
        expect_identical(output$provenance$core_version,result$provenance$core_version)
        rm(handle,leading);invisible(gc())
    }
})
