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
        #### Read VCF ####
        vcf <- VariantAnnotation::readVcf(file = path)
        #### Parse VCF ####
        sumstats_file <- vcf2df(v = vcf)  %>%
            dplyr::rename(CHROM = seqnames, 
                          BP = start) 
        #### Remove sample_id suffixes from columns ####
        # Might not be the best way to do this if there's multiple samples/GWAS
        sample_ids <- get_vcf_sample_ids(path = path)
        if(!is.null(sample_ids)){
            colnames(sumstats_file) <- gsub(
                paste(paste0("_",sample_ids,"$"),collapse = "|"),"",
                colnames(sumstats_file))
        } 
        #### Convert to data.tabel ####
        sumstats_file <- data.table::data.table(sumstats_file) 
        sumstats_file[,PARSED:=TRUE]
        sumstats_file
    })
    return(sumstats_file)
}
