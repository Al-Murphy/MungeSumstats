#' Ensure that CHR:BP aren't merged into 1 column
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
check_two_step_col <- function(sumstats_file, path){
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
  twoStepCol <- grep(".*:.*", row_of_data)
  #in case there are more than one column with ":", just take first one
  if (length(twoStepCol)>1){
    #sort to get most recent genome build by default (cols: SNP_hg19, SNP_hg18)
    keep_col <- sort(col_headers[twoStepCol],decreasing = TRUE)[1]
    drop_cols <- sort(col_headers[twoStepCol],decreasing = TRUE)[-1]
    msg <- paste0("Warning: Multiple columns in the sumstats file seem to ",
                  "relate to Chromosome:Base Pair position.\nThe column ",
                  keep_col," will be kept whereas the column(s) ",
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
    twoStepCol <- grep(".*:.*", row_of_data) # this should only have 1 col now
  }
  if (length(twoStepCol)) {
    # Convert the ':' into '\t'
    sumstats_file <-
      gsub(pattern = ":", replacement = "\t", x = sumstats_file)
    # Replace the column name with two names
    curColName <- col_headers[twoStepCol]
    # Write the new column headers to file
    new_first_line <-
      gsub(curColName,"CHR\tBP",paste(col_headers,collapse = "\t"))
    sumstats_file[1] <- new_first_line
    col_headers <- strsplit(new_first_line,"\t")[[1]]
    message(sprintf("Column %s has been replaced with CHR BP",curColName))
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
