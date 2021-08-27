#' Ensure valid tabular format
#'
#' @param header The summary statistics file for the GWAS
#' @return Whether the file is tabular
#' @keywords internal
check_tabular <- function(header) {
    ### If data.table can read in the header properly, it can read in the rest of the file.
    row_of_data <- colnames(data.table::fread(text = header))
    # check if readLines picked up headers as single char
    is_tabular <- length(row_of_data) != 1
    if (is_tabular) message("Tabular format detected.")
    return(is_tabular)
}
