a <- commandArgs(TRUE)
if(length(a)!=1L)stop('reference directory required')
m <- readRDS(file.path(a[[1]],'manifest.rds'))
paths <- c(fasta=m$fasta,stats::setNames(file.path(a[[1]],m$files),m$files))
out <- data.frame(kind='file',name=names(paths),path=unname(paths),
                  md5=unname(m$checksums[names(paths)]),
                  size=unname(m$sizes[names(paths)]),value='',stringsAsFactors=FALSE)
metadata <- data.frame(kind='metadata',name=c('build','version','core_version'),path='',md5='',size=0,
                       value=c(m$build,m$version,m$core_version),stringsAsFactors=FALSE)
write.table(rbind(metadata,out),stdout(),sep='\t',row.names=FALSE,quote=TRUE)
