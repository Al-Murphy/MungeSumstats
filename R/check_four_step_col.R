#' Ensure that CHR:BP:A2:A1 aren't merged into 1 column
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return The modified sumstats_file
check_four_step_col <- function(sumstats_file){
  #get col headers
  col_headers <- sumstats_file[1]
  col_headers <- strsplit(col_headers, "\t")[[1]]
  # Obtain a row of the actual data
  row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
  # Check if there is a column of data with CHR:BP:A2:A1 format
  fourStepCol <- grep(".*:.*:\\w:\\w",row_of_data)
  if(length(fourStepCol)){
    # Convert the ':' into '\t'
    sumstats_file <- 
      gsub(pattern = ":", replacement = "\t", x = sumstats_file)
    # Replace the column name with four names
    curColName <- col_headers[fourStepCol]
    # Write the new column headers to file
    first_line <- paste(col_headers, collapse = "\t")
    new_first_line <- 
      gsub(sprintf("^%s\\t|\\t%s\\t|\\t%s$",curColName,curColName,curColName),
           "CHR\tBP\tA2\tA1\t", paste(col_headers, collapse = "\t"))
    sumstats_file[1] <- new_first_line
    col_headers <- strsplit(new_first_line, "\t")[[1]]
    print(sprintf("Column %s has been replaced with CHR BP A2 A1", 
                  curColName))
    print(col_headers)
    
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}