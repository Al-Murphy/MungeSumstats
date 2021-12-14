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
    # Sometimes the N column is not all integers... so round it up
    col_headers <- names(sumstats_dt)
    if ("P" %in% col_headers) {
        # get smallest p-val - seems to change to character if < xe-300
        char_check <- FALSE
        num_check <- FALSE
        # Let negative p-values be handled by check_range_pval in the next step
        small_bool <- as.numeric(sumstats_dt$P) <= 5e-324 & 
                      as.numeric(sumstats_dt$P) > 0L
        small_n <-  sum(small_bool, na.rm = TRUE)
        
        if (is.numeric(sumstats_dt$P)) {
            if (small_n > 0) {
                num_check <- TRUE
            }
        } else { # char check
            max_minus_power <- max(
                as.numeric(gsub(".*-", "", sumstats_dt$P)),
                na.rm = TRUE)
            if (max_minus_power >= 324) {
                char_check <- TRUE
            }
        }   
        
        if (char_check | num_check) { 
            # check if any smaller or equal to 5e-324 limit
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
                sumstats_dt[small_bool, P := 0]
            } else {
                msg2 <- paste0(msg, "These will NOT be converted.")
                message(msg2)
            }
        }
    }
    return(list("sumstats_dt" = sumstats_dt))
}
