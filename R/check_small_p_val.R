#' Ensure that the non-negative p-values are not 5e-324 or lower, if so set to 0
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @inheritParams format_sumstats
#' @source 
#' \code{
#' sumstats_dt <- MungeSumstats:::formatted_example()
#' sumstats_dt$P[1:3] <- 5e-324
#' sumstats_dt$P[6:10] <- "5e-324"
#' sumstats <- check_small_p_val(sumstats_dt = sumstats_dt,
#'                               convert_small_p = TRUE, 
#'                               imputation_ind = TRUE)
#' }
#' @returns list containing sumstats_dt, 
#' the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table :=
check_small_p_val <- function(sumstats_dt,  
                              convert_small_p,
                              imputation_ind) {
    P <- convert_small_p_0 <- NULL
    col_headers <- names(sumstats_dt)
    if ("P" %in% col_headers) {
        p_num <- suppressWarnings(as.numeric(sumstats_dt$P))
        # Let negative p-values be handled by check_range_pval in the next step
        small_bool <- !is.na(p_num) & p_num <= 5e-324 & p_num > 0
        small_n <- sum(small_bool)
        # P read as text (fread keeps e.g. "1e-400" as character) can be below
        # the smallest number R can store (~4.9e-324), so as.numeric() gives 0
        underflow_bool <- rep(FALSE, length(p_num))
        if (!is.numeric(sumstats_dt$P)) {
            mantissa <- suppressWarnings(
                as.numeric(sub("[eE].*$", "", sumstats_dt$P)))
            underflow_bool <- !is.na(p_num) & p_num == 0 &
                !is.na(mantissa) & mantissa > 0
        }
        underflow_n <- sum(underflow_bool)

        if (underflow_n > 0) {
            message(
                formatC(underflow_n, big.mark = ","),
                " p-values are below ~4.9e-324, the smallest number R can ",
                "store, so they will be set to 0 whatever convert_small_p is. ",
                "Supply -log10 p-values to keep them."
            )
            if (imputation_ind) {
                sumstats_dt[underflow_bool, convert_small_p_0 := TRUE]
            }
        }
        if (small_n > 0) {
            msg <- paste0(
                formatC(small_n, big.mark = ","),
                " p-values are <=5e-324 which ",
                "LDSC/MAGMA may not be able to handle. "
            )
            if (convert_small_p) { # if user wants to correct
                msg2 <- paste0(msg, "These will be converted to 0.")
                message(msg2)
                # if users want edited snps, return information
                if (imputation_ind) {
                    sumstats_dt[small_bool, convert_small_p_0 := TRUE]
                }
                p_num[small_bool] <- 0
            } else {
                msg2 <- paste0(msg, "These will NOT be converted.")
                message(msg2)
            }
        }
        if (underflow_n > 0 || (small_n > 0 && convert_small_p)) {
            sumstats_dt[, P := p_num]
        }
    }
    return(list("sumstats_dt" = sumstats_dt))
}
