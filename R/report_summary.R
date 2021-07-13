#' Report info on current state of the summary statistics 
#' 
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @keywords internal  
report_summary <- function(sumstats_dt){
    message("Formatted summary statistics report:",
            "\n   - ", formatC(nrow(sumstats_dt), big.mark = ",")," rows",
            "\n   - ", formatC(length(unique(sumstats_dt$SNP)),big.mark = ",")," unique variants",
            "\n   - ", formatC(nrow(subset(sumstats_dt, P<5e-8)),big.mark = ",")," genome-wide significant variants (P<5e-8)",
            "\n   - ", formatC(length(unique(sumstats_dt$CHR)),big.mark = ",")," chromosomes"
            ) 
}