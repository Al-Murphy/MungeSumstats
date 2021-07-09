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
                           nThread=1){
    if(write_vcf){
        message("Writing in VCF format ==> ",check_save_out$save_path)
        # gr <- to_GRanges(sumstats_dt = sumstats_dt)
        vr <- to_VRanges(sumstats_dt = sumstats_dt)
        VariantAnnotation::writeVcf(vr, filename = check_save_out$save_path)
        
    } else {
        message("Writing in tabular format ==> ",check_save_out$save_path)
        data.table::fwrite(x = sumstats_dt,
                           file = check_save_out$save_path,
                           sep = check_save_out$sep, 
                           nThread = nThread) 
    } 
}