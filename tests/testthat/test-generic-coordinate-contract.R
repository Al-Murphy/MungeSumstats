test_that('complex generic ID inference requires source position even for one canonical event', {
    fixture <- .alt_orientation_reference(c('h95\t3\tA\tCG,ACCA,AACC\tevent95\told_event95',
                                            'unique\t1\tCCAG\tCCTG\tu\told_u'),
                                          c(h95='CCAACTCAGCTT',unique='CCAGG'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    for(chromosome in c('h95',NA_character_)) {
        x <- .alt_orientation_query(c('A','ACCA','T','TGGT'),c('ACCA','A','TGGT','T'),
                                    chr=chromosome,pos=NA_real_,id='old_event95')
        actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
        expect_identical(actual$source,x)
        .alt_orientation_expect_unresolved(actual$normalized,'source_coordinates_missing')
        expect_equal(actual$normalized$candidate_count,rep(3,4))
    }
    # One catalog event still admits two valid measurement anchors with opposite signs.
    for(chromosome in c('unique',NA_character_)) {
        x <- .alt_orientation_query(c('CCAG','CCTG','CTGG','CAGG'),
                                    c('CCTG','CCAG','CAGG','CTGG'),chr=chromosome,
                                    pos=NA_real_,id='old_u')
        actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
        expect_identical(actual$source,x)
        .alt_orientation_expect_unresolved(actual$normalized,'source_coordinates_missing')
        expect_equal(actual$normalized$candidate_count,rep(1,4))
    }
    absent <- .alt_orientation_query('A','AA',chr=NA_character_,pos=NA_real_,id='absent')
    result <- harmonise_sumstats(absent,fixture$reference)
    .alt_orientation_expect_unresolved(result$normalized,'source_coordinates_missing')
    expect_equal(result$normalized$candidate_count,0)
})

test_that('known source positions preserve distinct identities and opposite valid effect signs', {
    fixture <- .alt_orientation_reference(c('h95\t3\tA\tCG,ACCA,AACC\tevent95\t.',
                                            'unique\t1\tCCAG\tCCTG\tu\t.'),
                                          c(h95='CCAACTCAGCTT',unique='CCAGG'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    x <- .alt_orientation_query('A','ACCA',chr='h95',pos=c(3,4),id='event95')
    x$BETA <- .625;x$FRQ <- .137
    actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')$normalized
    expect_identical(actual$status,c('match','match'))
    expect_equal(actual$pos,c(1,3));expect_identical(actual$ref,c('C','A'))
    expect_identical(actual$alt,c('CCAC','AACC'));expect_equal(actual$effect,c(.625,.625))
    y <- .alt_orientation_query('CCAG','CCTG',chr='unique',pos=c(1,2),id='u')
    y$BETA <- .625;y$FRQ <- .137
    known <- harmonise_sumstats(y,fixture$reference,frequency_type='EAF')$normalized
    expect_identical(known$status,c('match','swapped_complement'))
    expect_equal(known$pos,c(3,3));expect_identical(known$ref,c('A','A'))
    expect_identical(known$alt,c('T','T'));expect_equal(known$effect,c(.625,-.625))
    expect_equal(known$eaf,c(.137,1-.137))
})

test_that('paired ID contexts survive coordinate candidate intersection and reference order', {
    rows <- c('pool\t2\tC\tT\tr\told_r','pool\t1\tAC\tGC,GT\tr\told_r')
    for(order in list(1:2,2:1)) {
        fixture <- .alt_orientation_reference(rows[order],c(pool='ACAAAAAAA'))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        for(id in c('r','old_r',NA_character_)) {
            x <- .alt_orientation_query('C','T',chr='pool',pos=2,id=id)
            actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
            expect_identical(actual$source,x)
            .alt_orientation_expect_unresolved(actual$normalized)
            expect_equal(actual$normalized$candidate_count,1)
            if(!is.na(id)) {
                x$CHR <- NA_character_;x$BP <- NA_real_
                unknown <- harmonise_sumstats(x,fixture$reference)$normalized
                .alt_orientation_expect_unresolved(unknown)
                expect_equal(unknown$candidate_count,3)
            }
        }
    }
})

test_that('closure keeps chromosome and paired REF context boundaries and trusted contracts', {
    fixture <- .alt_orientation_reference(c('pool\t2\tC\tT\tr\told_r',
                                            'foreign\t1\tAC\tGC,GT\tr\told_r'),
                                          c(pool='ACAAAAAAA',foreign='ACAAAAAAA'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    x <- .alt_orientation_query('C','T',chr='pool',pos=2,id='old_r')
    result <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')$normalized
    expect_identical(result$status,'match');expect_equal(result$effect,.3)
    expect_equal(result$candidate_count,1)
    columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',
                    source_ref='REF',source_alt='ALT',source_strand='STRAND')
    same <- .alt_orientation_reference(c('pool\t2\tC\tT\tr\told_r','pool\t1\tAC\tGC,GT\tr\told_r'),c(pool='ACAAAAAAA'))
    on.exit(unlink(same$root,recursive=TRUE),add=TRUE)
    x$REF <- 'C';x$ALT <- 'T';x$STRAND <- '+'
    declared <- harmonise_sumstats(x,same$reference,columns=columns,frequency_type='EAF')
    expect_identical(declared$source,x);expect_identical(declared$normalized$status,'match')
    expect_equal(declared$normalized$effect,.3)
    x$BP <- NA_real_
    .alt_orientation_expect_unresolved(harmonise_sumstats(x,same$reference,columns=columns)$normalized,'source_coordinates_missing')
})

test_that('single-base ID inference and unsupported complex marker behavior remain explicit', {
    fixture <- .alt_orientation_reference('1\t1\tA\tC\tr\told_r')
    on.exit(unlink(fixture$root,recursive=TRUE))
    snp <- .alt_orientation_query(c('A','C'),c('C','A'),chr=NA_character_,pos=NA_real_,id='old_r')
    actual <- harmonise_sumstats(snp,fixture$reference,frequency_type='EAF')$normalized
    expect_identical(actual$status,c('match','swap'))
    expect_equal(actual$effect,c(.3,-.3));expect_equal(actual$eaf,c(.2,.8))
    bad <- .alt_orientation_query(c('I','AA','N','<DEL>'),c('D','AA','A','A'),chr=NA_character_,pos=NA_real_,id='r')
    markers <- harmonise_sumstats(bad,fixture$reference)$normalized
    .alt_orientation_expect_unresolved(markers[1:3,],'unsupported_allele')
    .alt_orientation_expect_unresolved(markers[4,,drop=FALSE],'specialist_preserved')
})

test_that('newly required ID context lookups retain raw candidate cap failures', {
    fixture <- .alt_orientation_reference(c('pool\t2\tC\tT\tr\t.',
                                            'pool\t1\tAC\tGC,GT\tr\t.',
                                            'pool\t1\tN\tC,T\tother\t.',
                                            'pool\t1\tNN\tC,T\tother\t.',
                                            'pool\t1\tNNN\tC,T\tother\t.'),c(pool='ACAAAAAAA'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    expect_error(harmonise_sumstats(.alt_orientation_query('C','T',chr='pool',pos=2,id='r'),
                                   fixture$reference,candidate_limit=3),'candidate limit exceeded')
})

test_that('missing-position complex rows remain byte-preserved in streamed editions', {
    fixture <- .alt_orientation_reference('unique\t1\tCCAG\tCCTG\tu\t.',c(unique='CCAGG'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    x <- .alt_orientation_query(c('CCAG','CCTG'),c('CCTG','CCAG'),chr='unique',pos=NA_real_,id='u')
    path <- file.path(fixture$root,'input.tsv');write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
    output <- harmonise_file(path,fixture$reference,file.path(fixture$root,'edition'),frequency_type='EAF',batch_size=1)
    expect_equal(output$source_rows,2);expect_equal(output$normalized_rows,2)
    expect_identical(readBin(path,'raw',file.info(path)$size),readBin(file.path(output$path,'source-artifact'),'raw',file.info(path)$size))
    normalized <- read.delim(file.path(output$path,'normalized.tsv'),check.names=FALSE)
    .alt_orientation_expect_unresolved(normalized,'source_coordinates_missing')
    expect_equal(normalized$source_row,1:2)
})


test_that('partial historical aliases and exact-event coordinate bridges close current IDs', {
    fixture <- .alt_orientation_reference(c('pool\t2\tC\tT\tr\told_r',
                                            'pool\t1\tAC\tGC,GT\tr\t.'),c(pool='ACAAAAAAA'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    x <- .alt_orientation_query('C','T',chr='pool',pos=2,id='old_r')
    .alt_orientation_expect_unresolved(harmonise_sumstats(x,fixture$reference)$normalized)
    # ID-only r1 reaches r2 through the same normalized event, then r2's other
    # original context supplies the competing reverse T/G ALT-ALT contrast.
    bridge <- .alt_orientation_reference(c('pool\t2\tA\tC\tr1\t.',
                                           'pool\t2\tA\tC\tr2\t.',
                                           'pool\t1\tAA\tGT,GG\tr2\t.'),c(pool='AAAAAAAA'))
    on.exit(unlink(bridge$root,recursive=TRUE),add=TRUE)
    y <- .alt_orientation_query('A','C',chr=NA_character_,pos=NA_real_,id='r1')
    result <- harmonise_sumstats(y,bridge$reference)$normalized
    .alt_orientation_expect_unresolved(result);expect_equal(result$candidate_count,1)
})

test_that('closed-universe REF ALT hypotheses reject without promoting proof-only mappings', {
    rows <- c('pool\t1\tA\tG\tr\told_r',
              'pool\t1\tA\tG\tq\t.',
              'pool\t3\tA\tG\tq\t.')
    orders <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    for(order in orders) {
        fixture <- .alt_orientation_reference(rows[order],c(pool='AAAAA'))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        for(id in c('r','old_r')) {
            x <- .alt_orientation_query('A','G',chr=NA_character_,pos=NA_real_,id=id)
            unknown <- harmonise_sumstats(x,fixture$reference)$normalized
            .alt_orientation_expect_unresolved(unknown);expect_equal(unknown$candidate_count,1)
            x$CHR <- 'pool';x$BP <- 1
            known <- harmonise_sumstats(x,fixture$reference)$normalized
            expect_identical(known$status,'match');expect_equal(known$pos,1)
            expect_identical(known$id,'r');expect_equal(known$effect,.3)
            x$BP <- 3
            .alt_orientation_expect_unresolved(harmonise_sumstats(x,fixture$reference)$normalized,'coordinate_id_conflict')
        }
    }
    foreign <- .alt_orientation_reference(c('pool\t1\tA\tG\tr\t.',
                                            'pool\t1\tA\tG\tq\t.',
                                            'foreign\t3\tA\tG\tq\t.'),c(pool='AAAAA',foreign='AAAAA'))
    on.exit(unlink(foreign$root,recursive=TRUE),add=TRUE)
    x <- .alt_orientation_query('A','G',chr='pool',pos=NA_real_,id='r')
    expect_identical(harmonise_sumstats(x,foreign$reference)$normalized$status,'match')
})

test_that('current ID cycles terminate and fanout obeys explicit global closure caps', {
    cycle <- .alt_orientation_reference(c('pool\t1\tA\tC\tr1\tr2',
                                          'pool\t1\tA\tC\tr2\tr1'),c(pool='AAAAA'))
    on.exit(unlink(cycle$root,recursive=TRUE),add=TRUE)
    x <- .alt_orientation_query('A','C',chr='pool',pos=1,id='r1')
    actual <- harmonise_sumstats(x,cycle$reference,candidate_limit=2)$normalized
    expect_identical(actual$status,'match');expect_true(is.na(actual$id))
    expect_equal(actual$candidate_count,1);expect_equal(actual$effect,.3)
    fanout <- .alt_orientation_reference(c('pool\t1\tA\tC\tr1\told',
                                           'pool\t2\tA\tG\tr2\tr1',
                                           'pool\t3\tA\tC\tr3\tr2'),c(pool='AAAAA'))
    on.exit(unlink(fanout$root,recursive=TRUE),add=TRUE)
    x$SNP <- 'old'
    expect_error(harmonise_sumstats(x,fanout$reference,candidate_limit=2),'context closure limit exceeded')
    # The extra current ID on an unrelated known chromosome is filtered before
    # graph expansion, so it cannot exhaust the two-ID closure allowance.
    bounded <- .alt_orientation_reference(c('pool\t1\tA\tC\tr1\told',
                                            'pool\t2\tA\tG\tr2\tr1',
                                            'foreign\t3\tA\tC\tr3\tr2'),c(pool='AAAAA',foreign='AAAAA'))
    on.exit(unlink(bounded$root,recursive=TRUE),add=TRUE)
    expect_identical(harmonise_sumstats(x,bounded$reference,candidate_limit=2)$normalized$status,'match')
})
