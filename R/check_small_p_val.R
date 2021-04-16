#' Ensure that the p values are not 5e-324 or lower, if so set to 0
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param convert_small_p Binary, should p-values < 5e-324 be converted to 0? Small p-values pass the R limit and can cause errors with LDSC/MAGMA and should be converted. Default is TRUE.
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table :=
check_small_p_val <- function(sumstats_file, path, convert_small_p){
  P = NULL
  # Sometimes the N column is not all integers... so round it up
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if("P" %in% col_headers) {
    sumstats_dt <- data.table::fread(path)
    #get smallest p-val - seems to change to character if < xe-300
    char_check <- FALSE
    num_check <- FALSE
    if(is.numeric(sumstats_dt$P)){
      if(min(sumstats_dt$P)<=5e-324)
        num_check <- TRUE
    }
    else{#char check
      max_minus_power <- max(as.numeric(gsub(".*-","",sumstats_dt$P)))
      if(max_minus_power >= 324)
        char_check <- TRUE
    }
    if(char_check|num_check){ # check if any smaller or equal to 5e-324 limit
      msg <- paste0("There are existing p-values as low as 5e-324 which ",
                    "LDSC/MAGMA may not be able to handle. ")
      if (convert_small_p) { #if user wants to correct
        msg2 <- paste0(msg,"These will be converted to 0.")
        message(msg2)
        sumstats_dt[,P:=as.numeric(as.character(P))]
        data.table::fwrite(x=sumstats_dt, file=path, sep="\t")
        sumstats_file <- readLines(path)

        return(sumstats_file)
      }
      else{
        msg2 <-paste0(msg,"These will NOT be converted.")
        message(msg2)
      }
    }
  }
  return(sumstats_file)
}
