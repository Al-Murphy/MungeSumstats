#' Infer chrom format
#' 
#' Infer chromosome format of a VCF file.
#' @param path Path to VCF.
#' @param nrows Number of rows to sample.
#' @param verbose Print messages.
#' 
#' @keywords internal
#' @importFrom data.table fread
#' @returns logical: whether or not the chromosome contains "chr".
infer_chrom_format <- function(path,
                               nrows,
                               verbose=TRUE){
    messager("Inferring chr format.",v=verbose)
    header <- data.table::fread(text = readLines(con = path, 
                                                 n = nrows+1000),
                                skip = "#CHROM", 
                                nrows = nrows)
    has_chr <- grepl("chr",header[["#CHROM"]][1])
    return(has_chr)
}
