#' Read VCF: MarkerName column
#' 
#' Parse MarkerName/SNP column in VCF file.
#' 
#' @param sumstats_dt Summary stats data.table. 
#' @keywords internal
#' @returns Null output.
read_vcf_markername <- function(sumstats_dt){
    if(!"SNP" %in% colnames(sumstats_dt)){
        if ("MarkerName" %in% colnames(sumstats_dt)) { 
            messager("Renaming MarkerName as SNP.")
            data.table::setnames(x = sumstats_dt,
                                 old = "MarkerName", 
                                 new = "SNP")  
        } else if ("ID" %in% colnames(sumstats_dt) ) { 
            messager("Renaming ID as SNP.")
            data.table::setnames(x = sumstats_dt,
                                 old = "ID", 
                                 new = "SNP")  
        }
    } 
}