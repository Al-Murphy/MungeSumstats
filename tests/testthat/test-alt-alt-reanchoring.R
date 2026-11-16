test_that('common reference padding preserves unknown-strand contrast ambiguity', {
    forms <- list(c('1\t2\tA\tC,G,T\tshared\told_shared','GACAAAA'),
                  c('1\t1\tGA\tGC,GG,GT\tshared\told_shared','GACAAAA'),
                  c('1\t2\tAC\tCC,GC,TC\tshared\told_shared','GACAAAA'),
                  c('1\t1\tGAC\tGCC,GGC,GTC\tshared\told_shared','GACAAAA'))
    input <- .alt_orientation_query(c('A','G','T','C','A','C','T','G'),
                                    c('G','A','C','T','C','A','G','T'),pos=2)
    for(form in forms) {
        fixture <- .alt_orientation_reference(form[[1]],c(`1`=form[[2]]))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        for(mode in c('both','coordinate','id','alias')) {
            x <- input
            if(mode=='coordinate')x$SNP <- NA_character_
            if(mode=='id'){x$CHR <- NA_character_;x$BP <- NA_real_}
            if(mode=='alias')x$SNP <- 'old_shared'
            actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
            expect_identical(actual$source,x)
            .alt_orientation_expect_unresolved(actual$normalized)
        }
    }
})

test_that('jointly padded pure ALT contrasts retain their explicit uncertainty status', {
    forms <- c('1\t2\tA\tC,G\tshared\t.',
               '1\t1\tGA\tGC,GG\tshared\t.',
               '1\t2\tAC\tCC,GC\tshared\t.',
               '1\t1\tGAC\tGCC,GGC\tshared\t.')
    for(form in forms) {
        fixture <- .alt_orientation_reference(form,c(`1`='GACAAAA'))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        for(pos in c(2,NA_real_)) {
            x <- .alt_orientation_query('C','G',pos=pos)
            .alt_orientation_expect_unresolved(harmonise_sumstats(x,fixture$reference)$normalized,'alt_alt_contrast')
        }
    }
})

test_that('padded indel and MNV contexts retain all four orientation ambiguities', {
    forms <- list(
        list(row='1\t1\tA\tAT,T\tshared\t.',seq='AAAAAAA',pos=1,a=c('A','AT','T','AT'),b=c('AT','A','AT','T')),
        list(row='1\t1\tAA\tATA,TA\tshared\t.',seq='AAAAAAA',pos=1,a=c('A','AT','T','AT'),b=c('AT','A','AT','T')),
        list(row='1\t1\tCA\tCAT,CT\tshared\t.',seq='CAAAAAA',pos=2,a=c('A','AT','T','AT'),b=c('AT','A','AT','T')),
        list(row='1\t2\tAC\tAT,GT\tshared\t.',seq='GACAAAA',pos=2,a=c('AC','AT','GT','AT'),b=c('AT','AC','AT','GT')),
        list(row='1\t2\tACA\tATA,GTA\tshared\t.',seq='GACAAAA',pos=2,a=c('AC','AT','GT','AT'),b=c('AT','AC','AT','GT')),
        list(row='1\t1\tGAC\tGAT,GGT\tshared\t.',seq='GACAAAA',pos=2,a=c('AC','AT','GT','AT'),b=c('AT','AC','AT','GT')))
    for(form in forms) {
        fixture <- .alt_orientation_reference(form$row,c(`1`=form$seq))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        for(unknown in c(FALSE,TRUE)) {
            x <- .alt_orientation_query(form$a,form$b,pos=if(unknown)NA_real_ else form$pos)
            if(unknown)x$CHR <- NA_character_
            actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
            expect_identical(actual$source,x)
            if(unknown) {
                # Generic complex alleles require source coordinates; a site ID
                # alone cannot establish their measurement anchor or direction.
                .alt_orientation_expect_unresolved(actual$normalized,'source_coordinates_missing')
            } else .alt_orientation_expect_unresolved(actual$normalized)
        }
    }
})

test_that('shared alternate backgrounds remain part of the association contrast', {
    # Common G prefix is absent from REF. Trimming it must preserve the fact
    # that C/T was an ALT-ALT contrast, even for a numerically identical event.
    forms <- list(c('1\t1\tAA\tGC,GT,AG\tshared\t.','AAAAAAA'),
                  c('1\t1\tAC\tGC,GT,AT\tshared\t.','ACAAAAA'))
    for(form in forms) {
        fixture <- .alt_orientation_reference(form[[1]],c(`1`=form[[2]]))
        on.exit(unlink(fixture$root,recursive=TRUE),add=TRUE)
        for(pos in c(2,NA_real_)) {
            x <- .alt_orientation_query('C','T',pos=pos)
            .alt_orientation_expect_unresolved(harmonise_sumstats(x,fixture$reference)$normalized)
        }
    }
})

test_that('repeat-equivalent source anchors compare both ALT events in one REF span', {
    fixture <- .alt_orientation_reference('1\t1\tCAA\tCA,CAAA,CAAAA\tshared\t.',c(`1`='CAAAAAAAA'))
    on.exit(unlink(fixture$root,recursive=TRUE))
    x <- .alt_orientation_query('A','AAA',pos=3)
    actual <- harmonise_sumstats(x,fixture$reference,frequency_type='EAF')
    expect_identical(actual$source,x)
    .alt_orientation_expect_unresolved(actual$normalized)
})
