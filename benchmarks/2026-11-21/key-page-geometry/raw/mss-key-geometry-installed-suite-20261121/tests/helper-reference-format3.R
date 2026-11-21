.format3_fixture <- function() {
    root <- tempfile('typed-reference-');dir.create(root)
    fasta <- file.path(root,'reference.fa')
    writeLines(c('>1','CAAAAACGTACGT'),fasta,useBytes=TRUE)
    rows <- c('1\t8\tG\tA,T\trs雪\tancien_é;同義;old;old',
              '1\t4\tAA\tA\trsDel\told_del',
              '1\t8\tG\tA,T\trs雪\told;old',
              '1\t10\tA\tC,<DEL>\tthird\t.',
              '1\t10\tN\tG\tunsupported_only\told_bad',
              '1\t10\tA\t\t.\t.')
    variants <- file.path(root,'variants.tsv')
    writeLines(c('chr\tpos\tref\talt\tid\taliases',rows),variants,useBytes=TRUE)
    list(root=root,fasta=fasta,variants=variants,source_rows=length(rows),valid_alts=6L)
}
.format3_prepare <- function(fixture,format,name=paste0('ref-',format)) {
    prepare_reference(fixture$fasta,fixture$variants,file.path(fixture$root,name),
                      build='SYNTHETIC',version='typed-fixture',reference_format=paste0('mss-preservation-',format),batch_size=2L)
}
.format3_fields <- function(line) {
    fields <- strsplit(paste0(line,'\tEND'), '\t',fixed=TRUE)[[1]]
    fields[-length(fields)]
}
.format3_native <- function(path, fasta, cap=10000L) {
    mss_harmonise_cpp(path,fasta,'1',8,'r','G','A',.3,.2,'EAF',cap)
}
