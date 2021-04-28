#' Ensure that A1:A2 or A1/A2 or A1>A2 or A2>A1 aren't merged into 1 column
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
check_allele_merge <- function(sumstats_dt, path){
  #get col headers
  col_headers <- names(sumstats_dt)
  # Obtain a row of the actual data
  row_of_data <- as.character(sumstats_dt[1,])
  #criteria for alleles column are: contains spliter of : or / or > or <
  #and be three characters in length minus whitespaces: A/G
  twoAlleleCol <- grep(".*[/:><].*", gsub(" ", "", row_of_data, fixed = TRUE))
  corr_nchar <- (lapply(gsub(" ", "", row_of_data, fixed = TRUE), 
                          function(x) nchar(x))==3)[twoAlleleCol]
  #filter twoAlleleCol to those of nchar 3
  twoAlleleCol <-  twoAlleleCol[corr_nchar]
  #in case there are more than one column just take first one
  if (length(twoAlleleCol)>1 && sum(c("A1","A2") %in% col_headers)<=1){
    #sort to get most recent genome build (cols: alleles_hg19, alleles_hg18)
    keep_col <- sort(col_headers[twoAlleleCol],decreasing = TRUE)[1]
    drop_cols <- sort(col_headers[twoAlleleCol],decreasing = TRUE)[-1]
    msg <- paste0("Warning: Multiple columns in the sumstats file seem to ",
                  "relate to alleles A1>A2.\nThe column ",
                  keep_col," will be kept whereas the column(s) ",
                  drop_cols, " will be removed.\nIf this is not the correct ",
                  "column to keep, please remove all incorrect columns from ",
                  "those listed here before \nrunning `format_sumstats()`.")
    message(msg)
    #Get data without dropped
    sumstats_dt[, (drop_cols) := NULL]
    twoAlleleCol <- which(col_headers==keep_col)
  }
  if (length(twoAlleleCol)) {
    keep_col <- col_headers[twoAlleleCol]
    #now get character which is splitting the alleles in first row
    splitter <- 
      substring(sub(" ", "", row_of_data, fixed = TRUE)[twoAlleleCol],2,2)
    #split out col into separate values, keep names
    format <- strsplit(keep_col,":")[[1]]
    if(length(format)!=2)#check : and / and < and > in name
      format <- strsplit(keep_col,"/")[[1]]
    if(length(format)!=2)#check : and / and < and > in name
      format <- strsplit(keep_col,"<")[[1]]
    if(length(format)!=2)#check : and / and < and > in name
      format <- strsplit(keep_col,">")[[1]]
    if(length(format)!=2){#If neither found assign name
      format <- c("A1","A2")
      #lastly if symbol is < implies order is A2, A1 not A1,A2
      if(splitter=="<")
        format <- c("A2","A1")
    }
    sumstats_dt[, (format) := data.table::tstrsplit(get(keep_col),
                                                    split=splitter, fixed=TRUE)]
    #remove any leading/trailing whitespaces
    for(format_i in format)
      sumstats_dt[, (format_i) := gsub("\\s+$","", get(format_i), fixed = TRUE)]
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
