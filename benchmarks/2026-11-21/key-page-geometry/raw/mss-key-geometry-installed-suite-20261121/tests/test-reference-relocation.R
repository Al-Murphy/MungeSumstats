.relocation_fixture <- function(format = 2L, build = 'SYNTHETIC') {
    root <- tempfile('mss-transfer-');dir.create(root)
    original <- file.path(root, 'original reference');dir.create(original)
    fasta <- file.path(original, 'reference.fa')
    writeLines(c('>1', 'CAAAAACGTACGT'), fasta, useBytes=TRUE)
    fasta <- normalizePath(fasta,mustWork=TRUE)
    sites <- file.path(original, 'variants.tsv')
    writeLines(c('chr\tpos\tref\talt\tid\taliases',
                 '1\t8\tG\tA,T\trs雪\told;同義',
                 '1\t4\tAA\tA\tdel\told_del'), sites, useBytes=TRUE)
    ref <- prepare_reference(fasta,sites,file.path(original,'index'),build=build,version='transfer-fixture',
                             reference_format=paste0('mss-preservation-',format),contig_aliases=c(chr1='1'))
    input <- data.frame(CHR='chr1',BP=c('0008','8','4'),SNP=c('old','同義','old_del'),
                        A1=c('G','T','AA'),A2=c('T','G','A'),BETA=c('.300','-.300','.125'),
                        FRQ=c('.2','.8','.25'),untouched='raw:0001',stringsAsFactors=FALSE)
    baseline <- harmonise_sumstats(input,ref,frequency_type='EAF')
    files <- list.files(ref$path,full.names=TRUE)
    bundle <- stats::setNames(unname(tools::md5sum(files)),basename(files))
    list(root=root,original=original,fasta=fasta,reference=ref,input=input,baseline=baseline,bundle=bundle)
}
.relocation_move <- function(fixture) {
    # Streaming APIs can leave unreachable external pointers awaiting R's next
    # collection. Release those native streams before an offline transfer;
    # Windows correctly refuses moving a directory containing open files.
    invisible(gc())
    target <- file.path(fixture$root,'collaborator é copy')
    stopifnot(file.rename(fixture$original,target))
    list(path=file.path(target,'index'),fasta=file.path(target,'reference.fa'))
}
.relocation_bundle <- function(path) {
    files <- list.files(path,full.names=TRUE)
    stats::setNames(unname(tools::md5sum(files)),basename(files))
}

test_that('default open retains its existing scientific and receipt behavior', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        handle <- open_reference(f$reference)
        result <- harmonise_sumstats(f$input,handle,frequency_type='EAF')
        expect_identical(result,f$baseline)
        expect_null(handle$fasta_binding)
        expect_null(result$provenance$reference_fasta_binding)
        expect_identical(.relocation_bundle(f$reference$path),f$bundle)
        expect_equal(result$normalized$effect,c(.3,.3,.125))
        rm(handle);invisible(gc())
    }
})

test_that('explicit FASTA binds transferred immutable references without the original path', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        moved <- .relocation_move(f)
        expect_false(file.exists(f$fasta))
        expect_error(open_reference(moved$path),'checksum|size|exist')
        handle <- tryCatch(open_reference(moved$path,fasta=moved$fasta),error=function(e)e)
        expect_false(inherits(handle,'error'),info=if(inherits(handle,'error'))conditionMessage(handle) else '')
        if(inherits(handle,'error'))next
        expect_s3_class(handle,'mss_reference_handle')
        expect_identical(handle$fasta,normalizePath(moved$fasta,mustWork=TRUE))
        expect_identical(handle$fasta_binding$recorded_path,f$fasta)
        expect_identical(handle$fasta_binding$runtime_path,handle$fasta)
        expect_identical(handle$fasta_binding$implementation_version,'mss-fasta-binding-1')
        expect_identical(handle$fasta_binding$md5,f$reference$checksums[['fasta']])
        expect_identical(handle$fasta_binding$fai_md5,f$reference$checksums[['fasta.fai']])
        expect_identical(handle$checksums,f$reference$checksums)
        y <- harmonise_sumstats(f$input,handle,frequency_type='EAF')
        expect_identical(y$source,f$input)
        expect_identical(y$normalized,f$baseline$normalized)
        expect_identical(ncol(y$normalized),19L)
        expect_identical(y$provenance$implementation_version,'mss-preservation-7')
        expect_identical(y$provenance$core_version,paste0('mss-preservation-',format))
        expect_identical(y$provenance$reference_checksums,f$baseline$provenance$reference_checksums)
        expect_identical(y$provenance$reference_fasta_binding,handle$fasta_binding)
        expect_identical(.relocation_bundle(moved$path),f$bundle)
        expect_identical(readRDS(file.path(moved$path,'manifest.rds'))$fasta,f$fasta)
        rm(handle);invisible(gc())
    }
})

test_that('explicit FASTA never bypasses byte identity or path validation', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        expect_true('fasta' %in% names(formals(open_reference)))
        if(!'fasta' %in% names(formals(open_reference)))next
        wrong <- file.path(f$root,'same length wrong.fa')
        writeLines(c('>1','TAAAAACGTACGT'),wrong)
        expect_equal(file.info(wrong)$size,file.info(f$fasta)$size)
        expect_error(open_reference(f$reference,fasta=wrong),'checksum')
        for(path in list('',NA_character_,character(),c(f$fasta,f$fasta),file.path(f$root,'absent.fa'),f$root,1L))
            expect_error(open_reference(f$reference,fasta=path),'FASTA|fasta|path|file|exist')
        expect_identical(.relocation_bundle(f$reference$path),f$bundle)
    }
})

test_that('relocation verifies every unchanged index artifact including FAI and source', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        expect_true('fasta' %in% names(formals(open_reference)))
        if(!'fasta' %in% names(formals(open_reference)))next
        for(name in c('fasta.fai','records','source.tsv','keys','offsets')) {
            path <- file.path(f$root,paste0('corrupt-',name));dir.create(path)
            expect_true(all(file.copy(list.files(f$reference$path,full.names=TRUE),path)))
            file <- file.path(path,name);bytes <- readBin(file,'raw',file.info(file)$size)
            bytes[[1]] <- if(name=='fasta.fai') as.raw(10L) else as.raw(bitwXor(as.integer(bytes[[1]]),1L))
            writeBin(bytes,file)
            expect_error(open_reference(path,fasta=f$fasta),'checksum')
        }
        expect_identical(.relocation_bundle(f$reference$path),f$bundle)
    }
})

test_that('explicit rebinding rereads the manifest and owns a fresh physical FASTA handle', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        expect_true('fasta' %in% names(formals(open_reference)))
        if(!'fasta' %in% names(formals(open_reference)))next
        first <- open_reference(f$reference)
        second_path <- file.path(f$root,'second é physical copy.fa')
        expect_true(file.copy(f$fasta,second_path))
        poisoned <- first;poisoned$fasta <- 'wrong.runtime';poisoned$checksums[['fasta']] <- 'wrong.identity'
        second <- open_reference(poisoned,fasta=second_path)
        expect_false(identical(first$.native,second$.native))
        expect_identical(second$fasta,normalizePath(second_path,mustWork=TRUE))
        expect_identical(second$checksums,f$reference$checksums)
        expect_identical(harmonise_sumstats(f$input,second,frequency_type='EAF')$normalized,f$baseline$normalized)
        serialized <- file.path(f$root,'handle.rds');saveRDS(second,serialized)
        restored <- open_reference(readRDS(serialized),fasta=second_path)
        expect_identical(harmonise_sumstats(f$input,restored,frequency_type='EAF')$normalized,f$baseline$normalized)
        expect_false(identical(second$.native,restored$.native))
        expect_identical(second$fasta_binding,restored$fasta_binding)
        expect_identical(.relocation_bundle(f$reference$path),f$bundle)
        # Windows can deny timestamp writes while a native stream is open.
        # Assert the mutation outcome rather than treating a failed write as a
        # changed file. A successful timestamp change must invalidate reuse.
        before_stat <- file.info(second_path)[,c('size','mtime')]
        changed_time <- Sys.setFileTime(second_path,Sys.time()+10)
        after_stat <- file.info(second_path)[,c('size','mtime')]
        if (isTRUE(changed_time)) {
            expect_false(identical(before_stat,after_stat))
            expect_error(harmonise_sumstats(f$input,second,verify=FALSE),'checksum|stat')
        } else {
            expect_identical(after_stat,before_stat)
            expect_identical(harmonise_sumstats(f$input,second,frequency_type='EAF',verify=FALSE)$normalized,
                             f$baseline$normalized)
        }
        # source.tsv is not held by a native stream. A real size mutation tests
        # the complete binding's stat guard on every OS, including Windows.
        artifact_guard <- open_reference(f$reference,fasta=second_path)
        expect_identical(harmonise_sumstats(f$input,artifact_guard,frequency_type='EAF',verify=FALSE)$normalized,
                         f$baseline$normalized)
        artifact <- file.path(f$reference$path,'source.tsv')
        previous_size <- file.info(artifact)$size
        append_connection <- file(artifact,open='ab')
        writeBin(as.raw(10L),append_connection)
        close(append_connection)
        expect_equal(file.info(artifact)$size,previous_size+1)
        expect_error(harmonise_sumstats(f$input,artifact_guard,verify=FALSE),'checksum|stat')
        rm(first,second,restored,poisoned,artifact_guard);invisible(gc())
    }
})

test_that('streamed editions record physical FASTA binding without changing scientific fields', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        expect_true('fasta' %in% names(formals(open_reference)))
        if(!'fasta' %in% names(formals(open_reference)))next
        input <- file.path(f$root,'input.tsv');write.table(f$input,input,sep='\t',quote=TRUE,row.names=FALSE)
        base <- harmonise_file(input,f$reference,file.path(f$root,'before'),frequency_type='EAF')
        moved <- .relocation_move(f);handle <- open_reference(moved$path,fasta=moved$fasta)
        after <- harmonise_file(input,handle,file.path(f$root,'after'),frequency_type='EAF',batch_size=1L)
        for(name in c('normalized.tsv','source.tsv','source-artifact','reference_checksums.tsv'))
            expect_identical(unname(tools::md5sum(file.path(base$path,name))),unname(tools::md5sum(file.path(after$path,name))))
        expect_null(base$reference_fasta_binding)
        expect_identical(after$reference_fasta_binding,handle$fasta_binding)
        expect_identical(after$implementation_version,base$implementation_version)
        expect_identical(after$reference_checksums,base$reference_checksums)
        metadata <- read.delim(file.path(after$path,'metadata.tsv'),colClasses='character')
        value <- stats::setNames(metadata$value,metadata$key)
        expect_identical(value[['reference_fasta_recorded_path']],f$fasta)
        expect_identical(value[['reference_fasta_runtime_path']],handle$fasta)
        expect_identical(value[['reference_fasta_implementation_version']],'mss-fasta-binding-1')
        expect_identical(value[['reference_fasta_md5']],f$reference$checksums[['fasta']])
        expect_identical(value[['reference_fasta_fai_md5']],f$reference$checksums[['fasta.fai']])
        expect_identical(.relocation_bundle(moved$path),f$bundle)
        rm(handle);invisible(gc())
    }
})

test_that('reference-version refusals remain consistent with and without relocation', {
    for(format in 2:3) {
        f <- .relocation_fixture(format);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
        manifest <- file.path(f$reference$path,'manifest.rds');value <- readRDS(manifest)
        value$core_version <- 'mss-preservation-999';saveRDS(value,manifest)
        expect_error(open_reference(f$reference),'Unsupported reference schema')
        expect_true('fasta' %in% names(formals(open_reference)))
        if('fasta' %in% names(formals(open_reference)))
            expect_error(open_reference(f$reference,fasta=f$fasta),'Unsupported reference schema')
    }
})

test_that('projection receipts retain both explicit physical bindings and unchanged scientific identities', {
    for(format in 2:3) {
        s <- .relocation_fixture(format,'GRCh37');t <- .relocation_fixture(5L-format,'GRCh38')
        on.exit(unlink(s$root,recursive=TRUE),add=TRUE);on.exit(unlink(t$root,recursive=TRUE),add=TRUE)
        chain <- file.path(s$root,'identity.chain')
        writeLines(c('chain 100 1 13 + 0 13 1 13 + 0 13 1','13'),chain)
        projection <- prepare_projection(chain,s$reference,t$reference,'GRCh37','GRCh38')
        original_receipt <- projection$receipt
        before <- harmonise_sumstats(s$input,s$reference,frequency_type='EAF',projection=projection)
        input <- file.path(s$root,'input.tsv');write.table(s$input,input,sep='\t',quote=TRUE,row.names=FALSE)
        base <- harmonise_file(input,s$reference,file.path(s$root,'projection-before'),frequency_type='EAF',projection=projection)
        rm(projection);invisible(gc())
        ms <- .relocation_move(s);mt <- .relocation_move(t)
        hs <- open_reference(ms$path,fasta=ms$fasta);ht <- open_reference(mt$path,fasta=mt$fasta)
        projection <- prepare_projection(chain,hs,ht,'GRCh37','GRCh38')
        expect_identical(projection$source$fasta,normalizePath(ms$fasta,mustWork=TRUE))
        expect_identical(projection$target$fasta,normalizePath(mt$fasta,mustWork=TRUE))
        expect_identical(projection$receipt[names(original_receipt)],original_receipt)
        expect_identical(projection$receipt$source_fasta_binding,hs$fasta_binding)
        expect_identical(projection$receipt$target_fasta_binding,ht$fasta_binding)
        after <- harmonise_sumstats(s$input,hs,frequency_type='EAF',projection=projection)
        expect_identical(after$source,before$source)
        expect_identical(after$normalized,before$normalized)
        expect_identical(after$provenance$reference_checksums,before$provenance$reference_checksums)
        edition <- harmonise_file(input,hs,file.path(s$root,'projection-after'),frequency_type='EAF',projection=projection)
        for(name in c('normalized.tsv','source.tsv','source-artifact'))
            expect_identical(unname(tools::md5sum(file.path(base$path,name))),unname(tools::md5sum(file.path(edition$path,name))))
        receipt <- read.delim(file.path(edition$path,'projection.tsv'),colClasses='character')
        expect_identical(receipt$value[receipt$key=='target_fasta_binding.runtime_path'],ht$fasta)
        expect_identical(edition$projection_receipt$target_fasta_binding,ht$fasta_binding)
        expect_identical(.relocation_bundle(ms$path),s$bundle)
        expect_identical(.relocation_bundle(mt$path),t$bundle)
        rm(projection,hs,ht);invisible(gc())
    }
})
