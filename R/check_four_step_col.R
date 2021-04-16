#' Ensure that CHR:BP:A2:A1 aren't merged into 1 column
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
check_four_step_col <- function(sumstats_file, path){
  #get col headers
  col_headers <- sumstats_file[1]
  col_headers <- strsplit(col_headers, "\t")[[1]]
  # Obtain a row of the actual data
  row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
  # Check if there is a column of data with CHR:BP:A2:A1 format
  fourStepCol <- grep(".*:.*:\\w:\\w",row_of_data)
  #in case there are more than one column with ":", just take first one
  if (length(fourStepCol)>1){
    #sort to get most recent genome build by default (cols: SNP_hg19, SNP_hg18)
    keep_col <- sort(col_headers[fourStepCol],decreasing = TRUE)[1]
    drop_cols <- sort(col_headers[fourStepCol],decreasing = TRUE)[-1]
    msg <- paste0("Warning: Multiple columns in the sumstats file seem to ",
                  "relate to Chromosome:Base Pair position:A2:A1.\nThe column",
                  " ",keep_col," will be kept whereas the column(s) ",
                  drop_cols, " will be removed.\nIf this is not the correct ",
                  "column to keep, please remove all incorrect columns from ",
                  "those listed here before \nrunning `format_sumstats()`.")
    message(msg)
    #Get data without dropped
    write.table(x=utils::read.table(path)[-which(col_headers %in% drop_cols)],
                file=path, sep="\t", quote=FALSE, row.names = FALSE,
                col.names = FALSE)
    sumstats_file <- readLines(path)
    col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
    row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    fourStepCol <- grep(".*:.*:\\w:\\w",row_of_data)# should only have 1 col now
  }

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
           "\tCHR\tBP\tA2\tA1\t", paste(col_headers, collapse = "\t"))
    sumstats_file[1] <- new_first_line
    col_headers <- strsplit(new_first_line, "\t")[[1]]
    message(sprintf("Column %s has been replaced with CHR BP A2 A1",
                  curColName))

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
