#' Remove non-standard VCF columns
#'
#' @return sumstats_file
#'
#' @inheritParams format_sumstats
#' @inheritParams readLines 
#' @importFrom VariantAnnotation scanVcfHeader
#' @keywords internal 
remove_nonstandard_vcf_cols <- function(sample_id,
                                        sumstats_file,
                                        keep_extra_cols=FALSE){
    keep_cols <- c("CHROM","POS","ID", "REF","ALT","INFO", sample_id)
    colsToRemove <- names(sumstats_file)[!names(sumstats_file) %in% keep_cols]
    if(!is.null(colsToRemove) && keep_extra_cols==FALSE){
        message("Removing non-standard columns: ", paste(colsToRemove, collapse = ", "))
        sumstats_file[,(colsToRemove):=NULL]
    }
    return(sumstats_file)
}