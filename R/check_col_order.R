#' Ensure that the first three columns are SNP, CHR, BP in that order and
#' then A1, A2 if present
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics
#' data table object
#' @keywords internal
#' @importFrom data.table setcolorder
check_col_order <- function(sumstats_dt,
                            path) {
    col_headers <- names(sumstats_dt)
    # Use data tables for speed
    if (!sum(col_headers[seq_len(3)] == c("SNP", "CHR", "BP")) == 3 | (
        sum(c("A1", "A2") %in% col_headers) == 2 &
            !sum(col_headers[c(4, 5)] == c("A1", "A2")) == 2)) {
        msg <- paste0(
            "Reordering so first three column headers are SNP, CHR and",
            " BP in this order."
        )
        message(msg)
        whichSNP <- which(col_headers == "SNP")[1]
        whichCHR <- which(col_headers == "CHR")[1]
        whichBP <- which(col_headers == "BP")[1]
        # if A1 and A2 present, put them as four and five
        if (sum(c("A1", "A2") %in% col_headers) == 2) {
            message("Reordering so the fourth and fifth columns are A1 and A2.")
            whichA1 <- which(col_headers == "A1")[1]
            whichA2 <- which(col_headers == "A2")[1]
            otherCols <- setdiff(
                seq_len(length(col_headers)),
                c(whichSNP, whichCHR, whichBP, whichA1, whichA2)
            )
            data.table::setcolorder(sumstats_dt, c(
                whichSNP, whichCHR, whichBP, whichA1,
                whichA2, otherCols
            ))
        } else {
            otherCols <-
                setdiff(
                    seq_len(length(col_headers)),
                    c(whichSNP, whichCHR, whichBP)
                )
            data.table::setcolorder(
                sumstats_dt,
                c(whichSNP, whichCHR, whichBP, otherCols)
            )
        }
        return(list("sumstats_dt" = sumstats_dt))
    } else {
        return(list("sumstats_dt" = sumstats_dt))
    }
}
