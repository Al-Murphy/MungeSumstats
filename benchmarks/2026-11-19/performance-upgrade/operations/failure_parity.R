# Evaluate public installed APIs. Only intentional failures enter the comparison.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
.libPaths(c(args[[1]],args[[2]]), include.site=FALSE)
library(MungeSumstats,lib.loc=args[[1]])
scientific <- args[[3]]
root <- args[[4]]
stopifnot(!file.exists(root),dir.create(root))
check <- function(name, operation, pattern, edition) {
    err <- tryCatch({operation(); NULL}, error=identity)
    if (!inherits(err,"error") || !grepl(pattern,conditionMessage(err),ignore.case=TRUE) ||
        file.exists(edition)) stop("Invalid rejection or partial publication: ",name)
    # Neither scientific fixture paths nor their role-specific temporary
    # locations contribute to an error's meaning.
    message <- gsub(root,"<owned-root>",conditionMessage(err),fixed=TRUE)
    c(case=name,class=paste(class(err),collapse=","),error=message,published="FALSE")
}
records <- list()
for (format in 2:3) {
    refdir <- file.path(scientific,paste0("scientific_",format))
    fa <- file.path(scientific,"scientific.fa")
    handle <- open_reference(refdir)
    fixture <- function(name, raw, pattern, ...) {
        name <- paste0(name,"-f",format)
        input <- file.path(root,paste0(name,".tsv"))
        dest <- file.path(root,paste0(name,"-edition"))
        writeBin(raw,input)
        records[[name]] <<- check(name, function() harmonise_file(input,handle,dest,...),pattern,dest)
    }
    fixture("quote",charToRaw('CHR\tBP\nS\t"2'),"quote")
    fixture("width",charToRaw("CHR\tBP\nS\t2\textra\n"),"fields")
    fixture("limit",charToRaw("CHR\tBP\nS\t12345678901234567\n"),"max_record_bytes",max_record_bytes=12L)
    input <- charToRaw("CHR\tBP\tSNP\tA1\tA2\nD\t1\tdense\tA\tG\n")
    fixture("candidate",input,"candidate limit",candidate_limit=1L)
    fixture("receipt",input,"receipt",acquisition_receipt=list(rows=2L))
    archive <- file.path(root,paste0(format,".gz"))
    con <- gzfile(archive,"wb");writeBin(input,con);close(con)
    bytes <- readBin(archive,"raw",n=file.info(archive)$size)
    fixture("gzip",head(bytes,-8L),"gzip|compressed|read")
    bad <- file.path(root,paste0("corrupt-",format))
    stopifnot(dir.create(bad),all(file.copy(list.files(refdir,full.names=TRUE),bad)))
    con <- file(file.path(bad,"records"),"ab");writeBin(as.raw(0),con);close(con)
    dest <- file.path(root,"absent-edition")
    records[[paste0("corrupt",format)]] <- check(paste0("corrupt-f",format),
                 function() open_reference(bad,fasta=fa),"checksum|size",dest)
}
utils::write.table(do.call(rbind,records),file.path(root,"errors.tsv"),sep="\t",
                   row.names=FALSE,quote=TRUE,qmethod="double")
