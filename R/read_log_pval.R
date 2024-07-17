#' Read -log10 p-value column
#' 
#' Parse p-value column in VCF file.of other general -loq10 p-values
#' 
#' @param sumstats_dt Summary stats data.table. 
#' @param return_list Binary, whether to return the dt in a list or not - list 
#' is standard for the `format_sumstats()` function.
#' @inheritParams format_sumstats 
#' 
#' @returns Null output.
#' @keywords internal
#' @importFrom data.table setnames
read_log_pval <- function(sumstats_dt,mapping_file = sumstatsColHeaders,
                          return_list=TRUE){
    P <- NULL
    # Need to convert P-value, currently -log10
    #cgreate a check based on mapping file
    column_headers <- names(sumstats_dt) 
    column_headers_upper <- toupper(names(sumstats_dt))
    #get all p mappings 
    mapping_file_p <- mapping_file[mapping_file$Corrected=="P",]
    #now create all possible -log10 p-value
    log_p_col <- c(paste0("-LOG10_",mapping_file_p$Uncorrected),
                   paste0("LOG10_",mapping_file_p$Uncorrected),
                   paste0("-LOG10 ",mapping_file_p$Uncorrected),
                   paste0("LOG10 ",mapping_file_p$Uncorrected),
                   paste0("-LOG10-",mapping_file_p$Uncorrected),
                   paste0("LOG10-",mapping_file_p$Uncorrected),
                   paste0("-LOG10.",mapping_file_p$Uncorrected),
                   paste0("LOG10.",mapping_file_p$Uncorrected),
                   paste0("-LOG10",mapping_file_p$Uncorrected),
                   paste0("LOG10",mapping_file_p$Uncorrected),
                   paste0("LOG_",mapping_file_p$Uncorrected),
                   paste0("-L_",mapping_file_p$Uncorrected),
                   paste0("L_",mapping_file_p$Uncorrected),
                   paste0("LOG-",mapping_file_p$Uncorrected),
                   paste0("-L-",mapping_file_p$Uncorrected),
                   paste0("L-",mapping_file_p$Uncorrected),
                   paste0("LOG.",mapping_file_p$Uncorrected),
                   paste0("-L.",mapping_file_p$Uncorrected),
                   paste0("L.",mapping_file_p$Uncorrected),
                   paste0("-L",mapping_file_p$Uncorrected),
                   paste0("L",mapping_file_p$Uncorrected),
                   "LP")
    lp_found <- FALSE
    for(col_i in  column_headers_upper){
      if(isTRUE(any(col_i==log_p_col))){
        lp_found<-TRUE
        lp_col_up<-col_i
      }
    }
    if (isTRUE(lp_found)) {
        messager("sumstats has -log10 P-values; these will be",
                 "converted to unadjusted p-values in the 'P' column.")
        lp_col <- column_headers[[which(column_headers_upper %in% lp_col_up)]]
        #check if -log10 or log10
        if(max(sumstats_dt[,get(lp_col)])<=0){
          #log10 - less likely so more lenient check
          sumstats_dt[, P := 10^(as.numeric(get(lp_col)))]
        }else{
          #-log10 - more likely
          sumstats_dt[, P := 10^(-1 * as.numeric(get(lp_col)))]  
        }
        
    }
    #### Return format ####
    if(return_list){
      return(list("sumstats_dt" = sumstats_dt))
    }else {
      return(sumstats_dt)
    }
}