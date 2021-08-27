#' Check if save path and log folder is appropriate
#'
#' @return Corrected \code{save_path}, the file type, the separator, corrected
#' \code{log_folder},the log file extension
#'
#' @inheritParams format_sumstats
#' @keywords internal
check_save_path <- function(save_path, 
                            log_folder, 
                            log_folder_ind,
                            write_vcf = FALSE) {
    #### Add warning to users that temp files aren't actually saved ####
    if (dirname(save_path) == tempdir()) {
        message(
            "\n\n******::NOTE::******\n",
            "- Formatted results will be saved to `tempdir()` by default.\n",
            "- This means all formatted summary stats will be deleted",
            "upon ending the R session.\n",
            "- To keep formatted summary stats, change `save_path` ",
            "( e.g. `save_path=file.path('./formatted',basename(path))` ),",
            "  or make sure to copy files elsewhere after processing ",
            "( e.g. `file.copy(save_path, './formatted/' )`.\n",
            "********************",
            "\n"
        )
    }

    if (!is.logical(log_folder_ind)) {
        stop("log_folder_ind must be either TRUE or FALSE")
    }

    if (log_folder_ind && log_folder == tempdir()) {
        message(
            "******::NOTE::******\n",
            "- Log results will be saved to `tempdir()` by default.\n",
            "- This means all log data from the run will be ",
            "deleted upon ending the R session.\n",
            "- To keep it, change `log_folder` to an actual directory ",
            "(e.g. log_folder='./').\n",
            "********************",
            "\n"
        )
    }

    #### Do a bit of QC to get the full path ####
    ## Expand "~" into full path bc isn't recognized in certain envs(eg Python)
    save_path <- path.expand(save_path)
    log_folder <- path.expand(log_folder)
    ## Expand relative path "./" into absolute path bc it's less ambiguous
    save_path <- gsub("^[.]/", paste0(getwd(), "/"), save_path)
    log_folder <- gsub("^[.]/", paste0(getwd(), "/"), log_folder)

    suffixes <- supported_suffixes()
    if (is.null(save_path)) {
        save_path <- paste0(tempfile(), ".tsv.gz")
        # get extension of save path for log files
        extension <- ".tsv.gz"
        file_type <- "tempfile"
        sep <- "\t"
    } else {
        suffix_match <-
            vapply(suffixes, function(x) {
                grepl(x, tolower(save_path),
                    ignore.case = TRUE
                )
            },
            FUN.VALUE = logical(1)
            )
        if (sum(suffix_match) > 0) {
            file_type <- names(suffix_match)[suffix_match][1]
            # get extension of save path for log files
            extension <- file_type
            sep <- if (grepl(".csv", file_type)) "," else "\t"
        } else {
            if (write_vcf == FALSE) {
                stop(
                    "save_path file format not recognized: ", save_path,
                    "\nMust be one of: \n   ",
                    paste(suffixes, collapse = "\n   ")
                )
            }
        }
    }
    if (write_vcf) {
        save_path <- gsub(paste(suffixes, collapse = "|"), ".vcf.gz", save_path)
        sep <- "\t"
        file_type <- "vcf"
    } else {
        #### Account for mismatch between write_vcf and save_path ####
        suffixes.vcf <- supported_suffixes(
            tabular = FALSE,
            tabular_compressed = FALSE
        )
        if (any(endsWith(save_path, suffixes.vcf))) {
            message(
                "`write_vcf=FALSE` but `save_path` suggests VCF output ",
                "format. Switching output to tabular format (.tsv.gz)."
            )
            save_path <-
                gsub(paste(suffixes.vcf, collapse = "|"), ".tsv.gz", save_path)
            sep <- "\t"
            file_type <- ".tsv"
            # get extension of save path for log files
            # if output vcf, save log file .tsv.gz
            extension <- ".tsv.gz"
        }
        # get extension of save path for log files
        if (!exists("extension")) {
            # if output vcf, save log file type just .tsv.gz
            extension <- ".tsv.gz"
        }
    }
    #### Make sure dir exists
    if (is.character(save_path)) {
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    }
    dir.create(log_folder, showWarnings = FALSE, recursive = TRUE)

    message("Formatted summary statistics will be saved to ==> ", save_path)
    if (log_folder_ind) {
        message("Log data to be saved to ==> ", log_folder)
    }

    return(list(
        save_path = save_path,
        file_type = file_type,
        sep = sep,
        log_folder = log_folder,
        extension = extension
    ))
}
