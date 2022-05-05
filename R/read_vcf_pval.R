#' Read VCF: p-value column
#' 
#' Parse p-value column in VCF file.
#' 
#' @param sumstats_dt Summary stats data.table. 
#' 
#' @returns Null output.
#' @keywords internal
#' @importFrom data.table setnames
read_vcf_pval <- function(sumstats_dt){
    
    P <- LP <- NULL
    # Need to convert P-value, currently -log10
    if ("Pval" %in% colnames(sumstats_dt)) {
        data.table::setnames(sumstats_dt, "P","Pval")
    } else if ("LP" %in% names(sumstats_dt)) {
        messager("VCF file has -log10 P-values, these will be ",
                 "converted to unadjusted p-values in the 'P' column.")
        sumstats_dt[, P := 10^(-1 * as.numeric(LP))]
    } 
}