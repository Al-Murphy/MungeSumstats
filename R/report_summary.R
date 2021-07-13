#' Report info on current state of the summary statistics 
#' 
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @keywords internal  
report_summary <- function(sumstats_dt,
                           orig_dims=NULL){
    orig_dims_report <- if(!is.null(orig_dims)){
        paste0(" (",round(nrow(sumstats_dt)/orig_dims[1]*100,1),"% of original ",formatC(orig_dims[1],big.mark = ",", format = "fg")," rows)")
    } else {NULL};
    
    message("Summary statistics report:",
            "\n   - ", formatC(nrow(sumstats_dt), big.mark = ",", format = "fg")," rows",orig_dims_report,
            if("SNP" %in% colnames(sumstats_dt)){
                paste0("\n   - ", formatC(length(unique(sumstats_dt$SNP)),big.mark = ",", format = "fg")," unique variants")
            } else {NULL},
            if("P" %in% colnames(sumstats_dt)){
                paste0("\n   - ", formatC(nrow(subset(sumstats_dt, P<5e-8)),big.mark = ",", format = "fg")," genome-wide significant variants (P<5e-8)")
            } else {NULL},
            if("CHR" %in% colnames(sumstats_dt)){
                paste0("\n   - ", formatC(length(unique(sumstats_dt$CHR)),big.mark = ",")," chromosomes")
            } else {NULL}
    ) 
}