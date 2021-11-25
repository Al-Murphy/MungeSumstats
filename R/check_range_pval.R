#' Ensure that the p values are not >1 and if so set to 1
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @inheritParams format_sumstats
#' @source 
#' \code{
#' sumstats_dt <- MungeSumstats:::formatted_example()
#' sumstats_dt$P[1:3] <- 5
#' sumstats_dt$P[6:10] <- -5
#' sumstats <- check_range_p_val(sumstats_dt = sumstats_dt, 
#'                               convert_large_p = TRUE,
#'                               convert_neg_p = TRUE, 
#'                               imputation_ind = TRUE)

#' }
#' @returns list containing sumstats_dt, 
#' the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table :=
check_range_p_val <- function(sumstats_dt, 
                              convert_large_p,
                              convert_neg_p,
                              imputation_ind) {
    
    P <- convert_large_p_1 <- convert_neg_p_0 <- NULL
    # Sometimes the N column is not all integers... so round it up
    col_headers <- names(sumstats_dt)
    if ("P" %in% col_headers) {
        #### Ensure column is numeric ####
        sumstats_dt[,P:=as.numeric(P)] 
        
        #### Check for p-vals >1 ####
        large_bool <- sumstats_dt$P > 1
        large_n <- sum(large_bool, na.rm = TRUE)
        if (large_n > 0) {
            msg <- paste0(
                formatC(large_n, big.mark = ","),
                " p-values are >1 which ",
                "LDSC/MAGMA may not be able to handle. "
            ) 
            #### Convert p-values >1 to 1 ####
            if(convert_large_p){
                msg2 <- paste0(msg, "These will be converted to 1.")
                message(msg2) 
                # if users want edited snps, return information
                if (imputation_ind) {
                    sumstats_dt[large_bool, convert_large_p_1 := TRUE]
                }
                sumstats_dt[large_bool, P := 1]
            } else {
                msg2 <- paste0(msg, "These will NOT be converted.")
                message(msg2)
            } 
        }
        
        
        #### Check for p-vals <1 ####
        neg_bool <- sumstats_dt$P < 0
        neg_n <- sum(neg_bool, na.rm = TRUE)
        if (neg_n > 0) {
            msg3 <- paste0(
                formatC(neg_n, big.mark = ","),
                " p-values are <0 which ",
                "LDSC/MAGMA may not be able to handle. "
            ) 
            #### Convert p-values <0 to 0 ####
            if(convert_neg_p){
                msg4 <- paste0(msg3, "These will be converted to 0.")
                message(msg4) 
                # if users want edited snps, return information
                if (imputation_ind) {
                    sumstats_dt[neg_bool, convert_neg_p_0 := TRUE]
                }
                sumstats_dt[neg_bool, P := 0]
            } else {
                msg2 <- paste0(msg3, "These will NOT be converted.")
                message(msg2)
            } 
        }
    } 
    return(list("sumstats_dt" = sumstats_dt))
}
