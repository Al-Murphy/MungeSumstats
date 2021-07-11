#' Write sum stats file to disk
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @param check_save_out Output of \code{check_save_path()}.
#' @param write_vcf Whether to write as VCF (TRUE) or tabular file (FALSE).
#' @return \code{VRanges} object
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom VariantAnnotation makeVRangesFromGRanges
write_sumstats <- function(sumstats_dt,
                           check_save_out,
                           write_vcf=FALSE,
                           tabix_index=FALSE,
                           nThread=1){
    #### Make sure the directory actually exists
    dir.create(dirname(check_save_out$save_path), 
               showWarnings = FALSE, recursive = TRUE)
    
    if(write_vcf){  
        vr <- to_VRanges(sumstats_dt = sumstats_dt) 
        if(tabix_index) {
            suffixes <- supported_suffixes(tabular = FALSE, tabular_compressed = FALSE) 
            message("Writing in VCF format ==> ",gsub(paste(suffixes,collapse = "|"),
                                                      ".vcf.bgz",check_save_out$save_path))
            message("Compressing with bgzip and indexing with tabix.")
        } else { 
            message("Writing in VCF format ==> ",check_save_out$save_path) 
        }
        VariantAnnotation::writeVcf(vr, filename = check_save_out$save_path, 
                                    index = tabix_index)  
        
    } else {
        message("Writing in tabular format ==> ",check_save_out$save_path)
        data.table::fwrite(x = sumstats_dt,
                           file = check_save_out$save_path,
                           sep = check_save_out$sep, 
                           nThread = nThread) 
    } 
}