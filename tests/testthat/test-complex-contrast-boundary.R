test_that('a shared alternate background remains coherent after a repeat shift', {
    fixture <- .alt_orientation_reference('pool\t1\tCAA\tGA,GAAA,CAAAA\tr\told_r',c(pool='CAAAAAAAA'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    # At source position 4 the pair A/AAA can be REF/ALT insertion of two As.
    # The same observed strings also compare original ALT GA/GAAA on their
    # shared G background: deletion of one A versus insertion of one A, shifted
    # through the A repeat. That is an ALT-ALT association contrast.
    x <- .alt_orientation_query(c('A','AAA','T','TTT'),c('AAA','A','TTT','T'),chr='pool',pos=4,id='old_r')
    actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
    expect_identical(actual$source,x)
    .alt_orientation_expect_unresolved(actual$normalized,'source_allele_context_required')
})


test_that('single-base contrasts defer when an inserted alternate background cannot be excluded', {
    fixture <- .alt_orientation_reference('pool\t1\tAA\tGCA,GCG,AG\tr\t.',c(pool='AAAAAAA'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    x <- .alt_orientation_query(c('A','G','T','C'),c('G','A','C','T'),chr='pool',pos=2,id='r')
    actual <- harmonise_sumstats(x,fixture$reference)
    expect_identical(actual$source,x)
    .alt_orientation_expect_unresolved(actual$normalized,'source_allele_context_required')
})

test_that('length and base-pair exclusions retain ordinary supported contrasts', {
    indel <- .alt_orientation_reference('pool\t1\tA\tAT,ATTT\tr\t.',c(pool='AAAAAAA'))
    on.exit(unlink(indel$root,recursive=TRUE),add=TRUE)
    # ALT lengths differ by two; REF/ALT query lengths differ by one. They
    # cannot become this pair through shared trimming, padding or repeat shifts.
    x <- .alt_orientation_query('A','AT',chr='pool',pos=1,id='r')
    expect_identical(harmonise_sumstats(x,indel$reference)$normalized$status,'match')
    snp <- .alt_orientation_reference('pool\t1\tA\tC,G\tr\t.',c(pool='AAAAAAA'))
    on.exit(unlink(snp$root,recursive=TRUE),add=TRUE)
    x <- .alt_orientation_query(c('A','G'),c('G','A'),chr='pool',pos=1,id='r')
    expect_identical(harmonise_sumstats(x,snp$reference)$normalized$status,c('match','swap'))
    duplicate <- .alt_orientation_reference('pool\t1\tA\tC,C\tr\t.',c(pool='AAAAAAA'))
    on.exit(unlink(duplicate$root,recursive=TRUE),add=TRUE)
    x <- .alt_orientation_query('A','C',chr='pool',pos=1,id='r')
    expect_identical(harmonise_sumstats(x,duplicate$reference)$normalized$status,'match')
})

test_that('unknown alternate sequence defers generic effects and preserves explicit source contracts', {
    forms <- c('pool\t1\tA\tC,<DEL>\tr\t.',
               'pool\t1\tA\tC,.\tr\t.',
               'pool\t1\tA\tC,N\tr\t.')
    columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ',
                    source_ref='REF',source_alt='ALT',source_strand='STRAND')
    for(row in forms) {
        fixture <- .alt_orientation_reference(row,c(pool='AAAAAAA'))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        x <- .alt_orientation_query('A','C',chr='pool',pos=1,id='r')
        generic <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
        expect_identical(generic$source,x)
        .alt_orientation_expect_unresolved(generic$normalized,'source_allele_context_required')
        x$REF <- 'A';x$ALT <- 'C';x$STRAND <- '+'
        declared <- harmonise_sumstats(x,fixture$reference,columns=columns,frequency_type='EAF')$normalized
        expect_identical(declared$status,'match');expect_equal(declared$effect,.3)
    }
    composed <- .alt_orientation_reference('pool\t1\tCAA\tGA,GAAA,CAAAA\tr\t.',c(pool='CAAAAAAAA'))
    on.exit(unlink(composed$root,recursive=TRUE),add=TRUE)
    x <- .alt_orientation_query('A','AAA',chr='pool',pos=4,id='r')
    x$REF <- 'A';x$ALT <- 'AAA';x$STRAND <- '+'
    trusted <- harmonise_sumstats(x,composed$reference,columns=columns)$normalized
    expect_identical(trusted$status,'match');expect_equal(trusted$effect,.3)
})


test_that('jointly padded pure SNP contexts have one exact source anchor', {
    other_sites <- c('pool\t3\tA\tT,C\tr\t.',
                     'pool\t2\tAA\tAT,AC\tr\t.',
                     'pool\t3\tAA\tTA,CA\tr\t.',
                     'pool\t2\tAAA\tATA,ACA\tr\t.')
    for(other in other_sites)for(order in list(1:2,2:1)) {
        fixture <- .alt_orientation_reference(c('pool\t1\tA\tG\tr\t.',other)[order],c(pool='AAAAAAA'))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        x <- .alt_orientation_query(c('A','G','T','C'),c('G','A','C','T'),chr='pool',pos=1,id='r')
        actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
        expect_identical(actual$source,x)
        expect_identical(actual$normalized$status,c('match','swap','complement','swapped_complement'))
        expect_equal(actual$normalized$pos,rep(1,4))
        expect_equal(actual$normalized$effect,c(.3,-.3,.3,-.3))
        expect_equal(actual$normalized$eaf,c(.2,.8,.2,.8))
    }
})
