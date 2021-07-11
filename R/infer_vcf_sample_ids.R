#' Infer VCF sample ID(s)
#'
#' @inheritParams format_sumstats  
#' @return Tsample_id
#' @keywords internal 
infer_vcf_sample_ids <- function(sample_id=NULL,
                                 sumstats_file){
    if(is.null(sample_id)){
        idcol_index <- grep("FORMAT",colnames(sumstats_file),ignore.case = TRUE)
        if(any(length(colnames(sumstats_file))>=idcol_index)){
            sample_id <- colnames(sumstats_file)[seq(idcol_index[1]+1,length(colnames(sumstats_file)))] 
            message("sample ID(s) inferred from data colnames: ",paste(sample_id,collapse = ", "))
        } else {
            stop("Sample ID(s) could not be extracted from VCF header or inferred from data colnames.")
        }
    }
    return(sample_id)
}