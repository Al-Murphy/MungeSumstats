test_that('reference format2 remains the default and explicit format3 is opt-in', {
    f <- .format3_fixture();on.exit(unlink(f$root,recursive=TRUE))
    default <- prepare_reference(f$fasta,f$variants,file.path(f$root,'default'),build='SYNTHETIC',version='typed-fixture',batch_size=2L)
    expect_identical(default$core_version,'mss-preservation-2')
    explicit <- .format3_prepare(f,2)
    for(name in c('core-version','fasta.fai','keys','offsets','records','rejected.tsv','source.tsv')) {
        expect_identical(unname(tools::md5sum(file.path(default$path,name))),unname(tools::md5sum(file.path(explicit$path,name))))
    }
    typed <- .format3_prepare(f,3)
    expect_identical(typed$core_version,'mss-preservation-3')
    expect_error(prepare_reference(f$fasta,f$variants,file.path(f$root,'bad'),build='SYNTHETIC',version='v',reference_format='mss-preservation-4'),'reference_format|arg|format')
})

test_that('typed records factor contexts while retaining every raw occurrence and alias key', {
    f <- .format3_fixture();on.exit(unlink(f$root,recursive=TRUE))
    two <- .format3_prepare(f,2);three <- .format3_prepare(f,3)
    lines <- readLines(file.path(three$path,'records'),encoding='UTF-8')
    rows <- lapply(lines,.format3_fields);tags <- vapply(rows,`[[`,character(1),1)
    expect_equal(sum(tags=='S'),f$source_rows);expect_equal(sum(tags=='A'),f$valid_alts)
    expect_true(all(lengths(rows[tags=='A'])==7));expect_true(all(lengths(rows[tags=='S'])==5))
    expect_identical(rows[[which(tags=='S')[[1]]]][2:5],c('1','8','G','A,T'))
    keys <- lapply(c(two$path,three$path),function(path)sub('\t[^\t]+$','',readLines(file.path(path,'keys'),encoding='UTF-8')))
    expect_identical(table(keys[[1]]),table(keys[[2]]))
    expect_true(any(keys[[2]]=='I:同義'));expect_false(any(keys[[2]]=='I:old_bad'))
    expect_identical(readBin(f$variants,'raw',file.info(f$variants)$size),readBin(file.path(three$path,'source.tsv'),'raw',file.info(f$variants)$size))
})

test_that('typed and legacy references have identical public scientific rows', {
    f <- .format3_fixture();on.exit(unlink(f$root,recursive=TRUE))
    refs <- list(.format3_prepare(f,2),.format3_prepare(f,3))
    x <- data.frame(CHR='1',BP=c(8,8,8,4,10,10),SNP=c('old','同義','rs雪','rsDel','third','old_bad'),
                    A1=c('G','T','C','A','C','A'),A2=c('T','G','A','AA','A','C'),
                    BETA=c('.300','-.300','.300','-.300','.300','.300'),FRQ=c('.2','.8','.2','.8','.2','.2'),
                    untouched='keep:0001',stringsAsFactors=FALSE)
    values <- lapply(refs,function(reference)harmonise_sumstats(x,open_reference(reference),frequency_type='EAF',batch_size=1L))
    expect_identical(values[[1]]$source,x);expect_identical(values[[2]]$source,x)
    expect_identical(values[[1]]$normalized,values[[2]]$normalized)
    expect_identical(values[[2]]$provenance$implementation_version,'mss-preservation-7')
    expect_identical(values[[2]]$provenance$core_version,'mss-preservation-3')
    expect_identical(values[[2]]$provenance$reference_context_order,'semantic-s-context-v1')
})

test_that('typed file editions preserve original bytes and identify the data format', {
    f <- .format3_fixture();on.exit(unlink(f$root,recursive=TRUE))
    refs <- list(.format3_prepare(f,2),.format3_prepare(f,3))
    x <- data.frame(CHR='1',BP=c(8,4),SNP=c('old','old_del'),A1=c('G','AA'),A2=c('T','A'),BETA=c('.300','.125'),FRQ=c('.2','.25'))
    input <- file.path(f$root,'input.tsv');write.table(x,input,sep='\t',quote=FALSE,row.names=FALSE)
    outputs <- lapply(seq_along(refs),function(i)harmonise_file(input,refs[[i]],file.path(f$root,paste0('edition-',i)),frequency_type='EAF',batch_size=1L))
    expect_identical(readLines(file.path(outputs[[1]]$path,'normalized.tsv')),readLines(file.path(outputs[[2]]$path,'normalized.tsv')))
    for(output in outputs) {
        expect_equal(output$source_rows,nrow(x));expect_equal(output$normalized_rows,nrow(x))
        expect_identical(readBin(input,'raw',file.info(input)$size),readBin(file.path(output$path,'source-artifact'),'raw',file.info(input)$size))
    }
    metadata <- read.delim(file.path(outputs[[2]]$path,'metadata.tsv'),colClasses='character')
    expect_identical(metadata$value[metadata$key=='core_version'],'mss-preservation-3')
    expect_identical(metadata$value[metadata$key=='implementation_version'],'mss-preservation-7')
    expect_identical(metadata$value[metadata$key=='reference_context_order'],'semantic-s-context-v1')
})

test_that('all source target projection format combinations identify both formats honestly', {
    root <- tempfile('typed-projection-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
    fasta <- file.path(root,'reference.fa');writeLines(c('>1','ACGT'),fasta)
    variants <- file.path(root,'sites.tsv');writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t2\tC\tT\tr\t.'),variants)
    refs <- list()
    for(build in c('GRCh37','GRCh38'))for(format in 2:3) {
        refs[[paste(build,format)]] <- prepare_reference(fasta,variants,file.path(root,paste(build,format)),build=build,version='same-data',reference_format=paste0('mss-preservation-',format))
    }
    chain <- file.path(root,'identity.chain');writeLines(c('chain 100 1 4 + 0 4 1 4 + 0 4 1','4'),chain)
    x <- data.frame(CHR='1',BP=2,SNP='r',A1='C',A2='T',BETA=.3,FRQ=.2)
    expected <- NULL
    for(source_format in 2:3)for(target_format in 2:3) {
        p <- prepare_projection(chain,refs[[paste('GRCh37',source_format)]],refs[[paste('GRCh38',target_format)]],'GRCh37','GRCh38')
        y <- harmonise_sumstats(x,refs[[paste('GRCh37',source_format)]],frequency_type='EAF',projection=p)$normalized
        if(is.null(expected))expected <- y else expect_identical(y,expected)
        expect_identical(p$receipt$reference_format,if(source_format==target_format)paste0('mss-preservation-',source_format) else 'mixed')
        if(source_format==3||target_format==3) {
            expect_identical(p$receipt$source_reference_format,paste0('mss-preservation-',source_format))
            expect_identical(p$receipt$target_reference_format,paste0('mss-preservation-',target_format))
        }
    }
})

test_that('typed native records reject corruption without relying on manifest checksums', {
    root <- tempfile('typed-corruption-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
    fasta <- file.path(root,'reference.fa');writeLines(c('>1','CAAAAACGTACGT'),fasta)
    sites <- file.path(root,'sites.tsv');writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t8\tG\tA\tr\t.'),sites)
    ref <- prepare_reference(fasta,sites,file.path(root,'ref'),build='test',version='v',reference_format='mss-preservation-3')
    expect_identical(.format3_native(ref$path,fasta)$status,'match')
    original <- readBin(file.path(ref$path,'records'),'raw',file.info(file.path(ref$path,'records'))$size)
    original_text <- rawToChar(original)
    changes <- list(sub('A\t1\t8\tG\tA','X\t1\t8\tG\tA',original_text,fixed=TRUE),
                    sub('A\t1\t8\tG\tA','A\t1\t0\tG\tA',original_text,fixed=TRUE),
                    sub('A\t1\t8\tG\tA','A\t1\t8\tN\tA',original_text,fixed=TRUE),
                    sub('A\t1\t8\tG\tA','A\t1\t8\tG\tG',original_text,fixed=TRUE),
                    substr(original_text,1,nchar(original_text)-1L))
    for (i in seq_along(changes)) {
        dest <- file.path(root,paste0('bad-',i));dir.create(dest)
        expect_true(all(file.copy(list.files(ref$path,full.names=TRUE),dest)))
        writeBin(charToRaw(changes[[i]]),file.path(dest,'records'))
        expect_error(.format3_native(dest,fasta),'Typed reference')
    }
    for (byte in c(as.raw(0),as.raw(255))) {
        dest <- file.path(root,paste0('bad-byte-',as.integer(byte)));dir.create(dest)
        file.copy(list.files(ref$path,full.names=TRUE),dest)
        bytes <- original;bytes[[length(bytes)-2L]] <- byte
        writeBin(bytes,file.path(dest,'records'))
        expect_error(.format3_native(dest,fasta),'Typed reference')
    }
})

test_that('typed source framing and signed unsupported contexts are explicit', {
    f <- .format3_fixture();on.exit(unlink(f$root,recursive=TRUE))
    bytes <- readBin(f$variants,'raw',file.info(f$variants)$size)
    bytes[[length(bytes)-1L]] <- as.raw(255)
    bad <- file.path(f$root,'invalid-utf8.tsv');writeBin(bytes,bad)
    expect_error(prepare_reference(f$fasta,bad,file.path(f$root,'bad-utf8'),build='test',version='v',reference_format='mss-preservation-3'),'UTF-8')
    expect_false(file.exists(file.path(f$root,'bad-utf8')))
    signed <- file.path(f$root,'signed.tsv')
    writeLines(c('chr\tpos\tref\talt\tid\taliases','1\t0\tA\tG\tzero\t.','1\t-1\tA\tG\tnegative\t.'),signed)
    ref <- prepare_reference(f$fasta,signed,file.path(f$root,'signed'),build='test',version='v',reference_format='mss-preservation-3')
    expect_identical(readLines(file.path(ref$path,'records')),c('S\t1\t0\tA\tG','S\t1\t-1\tA\tG'))
    expect_identical(readLines(file.path(ref$path,'rejected.tsv')),c('row\tstatus','2\treference_mismatch','3\treference_mismatch'))
})
