#' Standardise the column headers in the Summary Statistics files
#'
#' Use a reference data table of common column header names (stored in
#' sumstatsColHeaders or user inputted mapping file) to convert them to a
#' standard set, i.e. chromosome -> CHR. This function does not check that all
#' the required column headers are present. The amended header is written
#' directly back into the file
#'
#' @inheritParams format_sumstats
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object
#' @keywords internal
#' @importFrom data.table setnames
standardise_sumstats_column_headers_crossplatform <-
    function(sumstats_dt,
             mapping_file = sumstatsColHeaders) {
        message("Standardising column headers.")
        message("First line of summary statistics file: ")
        msg <- paste0(names(sumstats_dt), split = "\t")
        message(msg)
        # first make all column headers upper case
        data.table::setnames(sumstats_dt, toupper(names(sumstats_dt)))
        column_headers <- names(sumstats_dt)
        # load synonym mapping - internal data no loading
        # Go through each and get correct spelling
        # allow for differing cases column names
        colnames(mapping_file) <- toupper(colnames(mapping_file))
        for (headerI in seq_len(nrow(mapping_file))) {
            un <- mapping_file[headerI, "UNCORRECTED"]
            cr <- mapping_file[headerI, "CORRECTED"]
            if (un %in% column_headers & (!cr %in% column_headers)) {
                data.table::setnames(sumstats_dt, un, cr)
            }
        }
        # Special case!! A0/A1 -> ref/alt so A1 flips meaning,
        # A0 is A* in mapping
        # but usually A1/A2 -> ref/alt so if A* found,
        # swap A1 to A2 and make A* -> A1
        new_headers <- colnames(sumstats_dt)
        if ("A*" %in% new_headers) {
            # if A1 and A2 also present need to rename A2
            if ("A1" %in% new_headers && "A2" %in% new_headers) {
                data.table::setnames(sumstats_dt, "A2", "A2_from_input")
            }
            # if A1 present change to A2, doesn't have to be, can be inputted
            data.table::setnames(sumstats_dt, "A1", "A2", skip_absent = TRUE)
            data.table::setnames(sumstats_dt, "A*", "A1")
        }
        return(list("sumstats_dt" = sumstats_dt))
    }
