#' Ensure that CHR:BP aren't merged into 1 column
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
check_two_step_col <- function(sumstats_dt, path){
  #get col headers
  col_headers <- names(sumstats_dt)
  # Obtain a row of the actual data
  row_of_data <- as.character(sumstats_dt[1,])
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
    sumstats_dt[, (drop_cols) := NULL]
    twoStepCol <- which(col_headers==keep_col)
  }
  if (length(twoStepCol)) {
    keep_col <- col_headers[twoStepCol]
    #split out col into separate values, keep names
    format <- strsplit(keep_col,":")[[1]]
    if(length(format)!=2)#check : and underscore in name
      format <- strsplit(keep_col,"_")[[1]]
    if(length(format)!=2)#If neither found assign name
      format <- c("CHR","BP")
    sumstats_dt[, (format) := data.table::tstrsplit(get(keep_col),
                                                    split=":", fixed=TRUE)]
    #remove combined column
    sumstats_dt[, (keep_col) := NULL]
    msg <- paste0("Column ",keep_col," has been separated into the columns ",
                  paste(format,collapse=", "))
    message(msg)
    
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
