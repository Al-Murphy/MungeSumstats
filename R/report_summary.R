#' Report info on current state of the summary statistics 
#' 
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @keywords internal  
report_summary <- function(sumstats_dt){
    message("Summary statistics report:\n",
            "\n   - ", nrow(sumstats_dt)," total variants",
            "\n   - ", nrow(subset(sumstats_dt, P<5e-8))," genome-wide significant variants (P<5e-8)",
            "\n   - ", nrow(sumstats_dt)," chromosomes"
            ) 
}