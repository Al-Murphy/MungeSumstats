read_vcf_data <- function(path,
                          nThread = 1,
                          tmpdir = tempdir(),
                          nrows  = Inf){
    start <- seqnames <- PARSED <- NULL;
    
    sumstats_file <- tryCatch(expr = {
        data.table::fread(
            input = path,
            nThread = nThread,
            sep = "\t",
            skip = "#CHR",
            tmpdir = tmpdir,
            nrows = nrows
        ) %>% dplyr::rename(CHROM = "#CHROM")
    }, error=function(e){
        vcf <- VariantAnnotation::readVcf(file = path)
        sumstats_file <- vcf2df(v = vcf)  %>%
            dplyr::rename(CHROM = seqnames, 
                          BP = start) 
        colnames(sumstats_file) <- gsub("_GWAS$","",colnames(sumstats_file))
        sumstats_file <- data.table::data.table(sumstats_file,
                                                keep.rownames = "ID") 
        sumstats_file[,PARSED:=TRUE]
        sumstats_file
    })
    return(sumstats_file)
}
