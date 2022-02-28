#' Simple function to ensure the new entry name to a list doesn't have the same
#' name as another entry
#'
#' @param name proposed name for the entry
#' @param log_files list of log file locations
#' @return a unique name (character)
#' @keywords internal
get_unique_name_log_file <- function(name, log_files) {
    used_names <- names(log_files)
    if (!name %in% used_names) { # if it isn't present return to be used
        return(name)
    }
    # else append number to end of name
    i <- 2
    need_unique_name <- TRUE
    while (need_unique_name) {
        proposed_name <- paste0(name, "_", i)
        if (!proposed_name %in% used_names) {
            name <- proposed_name
            need_unique_name <- FALSE
        }
        i <- i + 1
    }
    return(name)
}
