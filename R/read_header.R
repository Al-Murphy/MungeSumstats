#' Read in file header
#'
#' @inheritParams format_sumstats
#' @inheritParams readLines 
#' @importFrom VariantAnnotation scanVcfHeader
#' @keywords internal 
read_header <- function(path,
                        n=2, 
                        skip_vcf_metadata=FALSE){
    message("Reading header.")
    vcf_suffixes <- supported_suffixes(tabular = FALSE, tabular_compressed = FALSE) 
    if(startsWith(path,"https://gwas.mrcieu.ac.uk") | any(endsWith(path, vcf_suffixes))){
        if(skip_vcf_metadata){
            # gr <- GenomicRanges::GRanges(seqnames = 1, ranges = 1:10000)
            # param <- VariantAnnotation::ScanVcfParam(save_path, which=gr)
            # preview <- VariantAnnotation::readVcf(save_path, param = param)
            # preview <- VariantAnnotation::scanVcf(save_path, param=param)
            # data.table::fread(save_path, skip = "#", nrows = 10)
            header <- readLines(path, n = 100)
            i <- which(startsWith(header,"#CHR")) 
            header <- data.table::fread(text =  header[seq(i,i+n)] )   
        } else {
            header <- VariantAnnotation::scanVcfHeader(path) 
        }
    } else { 
        header <- readLines(con = path, n = n) 
        #Deal with strange, not recognised characters in header like '\' e.g 'xa6\xc2'
        header[[1]] <- iconv(enc2utf8(header[[1]]), sub="byte") 
    } 
    return(header)
}