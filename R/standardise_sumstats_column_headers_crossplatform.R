#' Standardise the column headers in the Summary Statistics files
#'
#' Use a reference data table of common column header names (stored in
#' \code{sumstatsColHeaders} or user inputted mapping file) to convert them to a
#' standard set, i.e. chromosome -> CHR. This function does not check that all
#' the required column headers are present. The amended header is written
#' directly back into the file
#'
#' @param uppercase_unmapped For columns that could not be identified in 
#' the \code{mapping_file}, return them in the same format they were input as
#'  (without forcing them to uppercase). 
#' @inheritParams format_sumstats 
#' @inheritParams compute_nsize
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object
#' @export
#' @importFrom data.table setnames rbindlist is.data.table
#' @examples
#' sumstats_dt <- data.table::fread(system.file("extdata", "eduAttainOkbay.txt",
#'                                              package = "MungeSumstats"))
#' sumstats_dt2 <- standardise_header(sumstats_dt=sumstats_dt)
standardise_header <- standardise_sumstats_column_headers_crossplatform <-
    function(sumstats_dt,
             mapping_file = sumstatsColHeaders,
             uppercase_unmapped=TRUE,
             return_list=TRUE) {
          
        stopifnot(
            "Mapping file must be a data.frame!" =
            !is.data.table(sumstatsColHeaders)
        )
        message("Standardising column headers.")
        message("First line of summary statistics file: ")
        msg <- paste0(names(sumstats_dt), split = "\t")
        message(msg)
        #### Store original colnames #####
        ## In case we want to use the original casing 
        ## IMPORTANT! Must use copy() function or else this vector will get 
        ## changed when the data.table it comes from gets changed.
        names_og <- data.table::copy(names(sumstats_dt))
        #### first make all column headers upper case ####
        data.table::setnames(sumstats_dt, toupper(names(sumstats_dt))) 
        column_headers <- names(sumstats_dt) 
        # load synonym mapping - internal data no loading
        # Go through each and get correct spelling
        # allow for differing cases column names
        colnames(mapping_file) <- toupper(colnames(mapping_file))
        #get names for allele mared eff/frq columns
        eff_frq_allele_matches <- get_eff_frq_allele_combns()
        mapping_file <- rbind(mapping_file,
                              as.data.frame(eff_frq_allele_matches))
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
        #### Return mapped columns to original casing if requested ####
        if(!uppercase_unmapped){
            unmapped <- colnames(sumstats_dt)[
                !colnames(sumstats_dt) %in% mapping_file$CORRECTED
            ]
            unmapped_og <- names_og[toupper(names_og) %in% unmapped]
            if(length(unmapped_og)>0){
                msg <- paste("Returning unmapped column names",
                             "without making them uppercase.")
                message(msg)
                data.table::setnames(sumstats_dt, unmapped, unmapped_og)
            } 
        }
        #### Return format ####
        if(return_list){
            return(list("sumstats_dt" = sumstats_dt))
        }else {
            return(sumstats_dt)
        } 
    }
