#' Ensure that CHR:BP aren't merged into 1 column
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return The modified sumstats_file
check_two_step_col <- function(sumstats_file){
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
  twoStepCol <- grep(".*:.*", row_of_data)
  if (length(twoStepCol)) {
    # Convert the ':' into '\t'
    sumstats_file <- 
      gsub(pattern = ":", replacement = "\t", x = sumstats_file)
    # Replace the column name with four names
    curColName <- col_headers[twoStepCol]
    # Write the new column headers to file
    new_first_line <- 
      gsub(curColName,"CHR\tBP",paste(col_headers,collapse = "\t"))
    sumstats_file[1] <- new_first_line
    col_headers <- strsplit(new_first_line,"\t")[[1]]
    print(sprintf("Column %s has been replaced with CHR BP",curColName))
    print(col_headers)
    
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}  