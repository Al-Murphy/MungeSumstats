#' Ensure that the p values are not 5e-324 or lower, if so set to 0
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table :=
check_small_p_val <- function(sumstats_file, path){
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
      max_minus_power <- max(gsub(".*-","",sumstats_dt$P))
      if(max_minus_power >= 324)
        char_check <- TRUE
    }
    if(char_check|num_check){ # check if any smaller or equal to 5e-324 limit
      msg <- paste0("MAGMA may not handle P-values as low as 5e-324. Do you ",
                    "want MAGMA.celltyping to \nconvert any existing ones",
                    "to zeroes? Type 0 for NO, 1 for YES: ")
      choice <- 2
      while(!choice %in% c(0,1)){
        choice <- readline(msg)
        if(!choice %in% c(0,1)){
          message(paste0(choice," is not a valid option. Please try again"))
        }
      }
      if (as.logical(as.numeric(readline(msg)))) { #if user wants to correct
        sumstats_dt[,P:=as.numeric(as.character(P))]
        data.table::fwrite(x=sumstats_dt, file=path, sep="\t")
        sumstats_file <- readLines(path)

        return(sumstats_file)
      }
    }
  }
  return(sumstats_file)
}
