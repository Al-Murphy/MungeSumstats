#' Ensure that CHR:BP:A2:A1 aren't merged into 1 column
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified
#' summary statistics data table object
#' @keywords internal
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
check_four_step_col <- function(sumstats_dt, path) {
    A2 <- A1 <- NULL
    # get col headers
    col_headers <- names(sumstats_dt)
    # Obtain a row of the actual data
    row_of_data <- as.character(sumstats_dt[1, ])
    # Check if there is a column of data with CHR:BP:A2:A1 format
    fourStepCol <- grep(".*:.*:\\w:\\w", row_of_data)
    # in case there are more than one column with ":", just take first one
    if (length(fourStepCol) > 1) {
        # sort to get most recent genome build by default
        # (cols: SNP_hg19, SNP_hg18)
        keep_col <- sort(col_headers[fourStepCol], decreasing = TRUE)[1]
        drop_cols <- sort(col_headers[fourStepCol], decreasing = TRUE)[-1]
        msg <- paste0(
            "Warning: Multiple columns in the sumstats file seem to ",
            "relate to Chromosome:Base Pair position:A2:A1.\nThe column",
            " ", keep_col, " will be kept whereas the column(s) ",
            drop_cols, " will be removed.\nIf this is not the correct ",
            "column to keep, please remove all incorrect columns from ",
            "those listed here before \nrunning `format_sumstats()`."
        )
        message(msg)
        # Get data without dropped
        sumstats_dt[, (drop_cols) := NULL]
        fourStepCol <- which(col_headers == keep_col)
    }
    if (length(fourStepCol)) {
        keep_col <- col_headers[fourStepCol]
        # split out col into separate values, keep names
        format <- strsplit(keep_col, ":")[[1]]
        if (length(format) != 4) { # check : and underscore in name
            format <- strsplit(keep_col, "_")[[1]]
        }
        if (length(format) != 4) { # If neither found
            # first check if allele col exists and assign based on that
            fourStepCol_val <- row_of_data[fourStepCol]
            if (sum("A2" %in% col_headers) == 1 && 
                sum("A1" %in% col_headers) == 1){
              A2_val <- sumstats_dt[1,A2]
              A1_val <- sumstats_dt[1,A1]
              split_fourStepCol <- strsplit(fourStepCol_val,split=":")[[1]]
              A1_fnd <- any(split_fourStepCol==A1_val)
              A2_fnd <- any(split_fourStepCol==A2_val)
              A1_ind <- which(split_fourStepCol==A1_val)
              A2_ind <- which(split_fourStepCol==A2_val)
            }
            #make sure you got a hit for A1 and A2 vals
            if(sum("A2" %in% col_headers) == 1 && 
               sum("A1" %in% col_headers) == 1 && 
               A1_fnd && A2_fnd){
              #bought A1 and A2 are present and we know their position
              format <- vector(mode="character", length=4)
              othrs <- c('CHR','BP')
              j <- 1
              for (i in seq_len(4)){
                if (i == A1_ind){
                  format[[i]] <- 'A1'
                } else if (i == A2_ind){
                  format[[i]] <- 'A2'
                } else{
                  format[[i]] <- othrs[[j]]
                  j = j+1
                }
              }
            } else if(sum("A2" %in% col_headers) == 1){
              A2_val <- sumstats_dt[1,A2]
              split_fourStepCol <- strsplit(fourStepCol_val,split=":")[[1]]
              A2_fnd <- any(split_fourStepCol==A2_val)
              A2_ind <- which(split_fourStepCol==A2_val)
              if (A2_ind == 4){
                format <- c("CHR", "BP", "A1", "A2")
              } else if(A2_ind == 3){
                format <- c("CHR", "BP", "A2", "A1")
              } else if(A2_ind == 2){
                format <- c("CHR", "A2","BP", "A1")
              } else if(A2_ind == 1){
                format <- c("A2", "CHR","BP", "A1")
              } else{
                #not found....
                format <- c("CHR", "BP", "A2", "A1")
              }
            } else if(sum("A1" %in% col_headers) == 1){
              A1_val <- sumstats_dt[1,A1]
              split_fourStepCol <- strsplit(fourStepCol_val,split=":")[[1]]
              A1_fnd <- any(split_fourStepCol==A1_val)
              A1_ind <- which(split_fourStepCol==A1_val)
              if (A1_ind == 4){
                format <- c("CHR", "BP", "A2", "A1")
              } else if(A1_ind == 3){
                format <- c("CHR", "BP", "A1", "A2")
              } else if(A1_ind == 2){
                format <- c("CHR", "A1","BP", "A2")
              } else if(A1_ind == 1){
                format <- c("A1", "CHR","BP", "A2")
              } else{
                #not found....
                format <- c("CHR", "BP", "A2", "A1")
              }
            }
            else{
              #otherwise, just use this as default order
              format <- c("CHR", "BP", "A1", "A2")
            }
            
        }
        sumstats_dt[, (format) := data.table::tstrsplit(get(keep_col),
            split = ":", fixed = TRUE,
            type.convert = TRUE
        )]
        # remove combined column
        sumstats_dt[, (keep_col) := NULL]
        msg <- paste0(
            "Column ", keep_col, " has been separated into the columns ",
            paste(format, collapse = ", "),"\nIf this is the incorrect ",
            "format for the column, update the column name to the correct ",
            "format e.g.`CHR:BP:A2:A1` and format_sumstats()."
        )
        message(msg)

        return(list("sumstats_dt" = sumstats_dt))
    } else {
        return(list("sumstats_dt" = sumstats_dt))
    }
}
