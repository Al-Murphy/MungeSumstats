#' Read in file header
#'
#' @return First \code{n} lines of the VCF header
#'
#' @inheritParams format_sumstats
#' @param n integer. The (maximal) number of lines to read. Negative values
#' indicate that one should read up to the end of input on the connection.
#' @param skip_vcf_metadata logical, should VCF metadata be ignored
# @inheritParams readLines
#' @importFrom VariantAnnotation scanVcfHeader
#' @keywords internal
read_header <- function(path,
                        n = 2,
                        skip_vcf_metadata = FALSE) {
    message("Reading header.")
    vcf_suffixes <- supported_suffixes(tabular = FALSE,
                                       tabular_compressed = FALSE)
    if (startsWith(path, "https://gwas.mrcieu.ac.uk") | 
        any(endsWith(path, vcf_suffixes))) {
        #### Read VCFs ####
        if (skip_vcf_metadata) {
            # gr <- GenomicRanges::GRanges(seqnames = 1, ranges = 1:10000)
            # param <- VariantAnnotation::ScanVcfParam(save_path, which=gr)
            # preview <- VariantAnnotation::readVcf(save_path, param = param)
            # preview <- VariantAnnotation::scanVcf(save_path, param=param)
            # data.table::fread(save_path, skip = "#", nrows = 10)
            header <- readLines(path, n = 100)
            i <- which(startsWith(header, "#CHR"))
            header <- data.table::fread(text = header[seq(i, i + n)],
                                        nThread = 1)
        } else {
            error_return <-
                tryCatch(header <-VariantAnnotation::scanVcfHeader(path),
                error = function(e) e,
                warning = function(w) w
                )
            if(is(error_return, "error")){
                #variant annotation can fail if so use brut force
                header <- readLines(path, n = 100)
                i <- which(startsWith(header, "#CHR"))
                header <- data.table::fread(text = header[seq(i, i + n)],
                                        nThread = 1)
            }    
        }
    } else if (endsWith(path,".bgz")){
        #### Read tabix-indexed tabular #### 
        # data.table::fread currently can't handle bgzipped files by default
        # 
        # header <- seqminer::tabix.read.header(tabixFile = path)$header 
        # header <- Rsamtools::headerTabix(file = path)$header
        # header <- rep(header,n)
        # header <- colnames(data.table::fread(text = header)) 
        header <- data.table::fread(cmd = paste("gunzip -c ",path),
                                    nrows = n) 
    } else if(startsWith(path, "https://")){
        #### Read generic remote sumstats ####
        header <- data.table::fread(path,nrows = n) 
    } else {
        #### Read tabular #### 
        header <- readLines(con = path,
                            n = n)
        # Deal with strange, not recognised characters in header
        # like '\' e.g 'xa6\xc2'
        header[[1]] <- iconv(enc2utf8(header[[1]]), sub = "byte")
        if(length(header)<2) header <- rep(header,2)
        header <- data.table::fread(text = header)
    }
    return(header)
}
