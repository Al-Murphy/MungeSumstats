#' Read VCF: INFO column
#' 
#' Parse INFO column in VCF file.
#' 
#' @param sumstats_dt Summary stats data.table. 
#' @keywords internal
#' @returns Null output.
#' @importFrom data.table setnames
read_vcf_info <- function(sumstats_dt){
    INFO <- NULL;
    
    if((!"INFO" %in% names(sumstats_dt)) && 
       ("SI" %in% names(sumstats_dt))){
        messager("Renaming SI column as INFO.")
        data.table::setnames(sumstats_dt, "SI","INFO")
    } 
    if("INFO" %in% names(sumstats_dt)){
        messager("Formatting INFO column.")
        sumstats_dt[INFO == ".", INFO := 0]
        sumstats_dt[, INFO := as.numeric(INFO)]
    } else {
        messager("No INFO (SI) column detected.")
    }
}
