#' Ensure that there is at least one signed column in summary statistics file
#'
#' @param sumstats_dt data table obj of the summary statistics 
#' file for the GWAS
#' @return null
#' @keywords internal
# check_signed_col <- function(sumstats_dt) {
#     col_headers <- names(sumstats_dt)
#     signed_stat_column_names <- c("Z", "OR", "BETA",
#                                   "LOG_ODDS", "SIGNED_SUMSTAT")
#     stp_msg <- paste0(
#         "ERROR: cannot find a column name representing signed ",
#         "statistic in GWAS sumstats file:\n",
#         "'Z','OR','BETA','LOG_ODDS','SIGNED_SUMSTAT'"
#     )
#     if (sum(signed_stat_column_names %in% col_headers) < 1) {
#         stop(stp_msg)
#     }
# }

check_signed_col <- function(sumstats_dt) {
    col_headers <- names(sumstats_dt)
    signed_stat_column_names <- c("Z", "OR", "BETA",
                                  "LOG_ODDS", "SIGNED_SUMSTAT")
    
    stp_msg <- paste0(
        "ERROR: cannot find a column name representing signed ",
        "statistic in GWAS sumstats file:\n",
        "'Z','OR','BETA','LOG_ODDS','SIGNED_SUMSTAT'"
    )
    
    msg <- "The sumstats Beta column is not present... "
    
    if ("BETA" %nin% col_headers) {
        
        if ("BETA" %nin% col_headers & "OR" %in% col_headers) {
            message(paste0(msg,"Deriving BETA from OR"))
            sumstats_dt[,BETA := log(OR)]
        } else if ("Z" %in% col_headers & "SE" %in% col_headers) {
            message(paste0(msg,"Deriving BETA from Z and SE"))
            sumstats_dt[,BETA := Z * SE]
        } else if ("Z" %in% col_headers & "N" %in% col_headers & "FRQ" %in% col_headers){
            message(paste0(msg,"Deriving BETA from Z, N, and FRQ"))
            sumstats_dt[,BETA := Z/sqrt(2*FRQ*(1-FRQ)*(N+Z^2))]
        } else if ("Z" %in% col_headers & "N" %in% col_headers & "P" %in% col_headers){
            message(paste0(msg,"Deriving BETA from Z, N, and P"))
            stop("Not yet implemented")
            # sumstats_dt[,BETA := Z/sqrt(chdtri(N, P))]
        } else if (sum(signed_stat_column_names %in% col_headers) < 1) {
            stop(stp_msg)
        }
        return(list("sumstats_dt" = sumstats_dt))
    }
    return(list("sumstats_dt" = sumstats_dt))
}
