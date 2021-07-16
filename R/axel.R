
#' axel downloader
#'
#' R wrapper for axel, which enables multi-threaded download of a single large file.
#' 
#' @return Path where the file has been downloaded
#' 
#' @inheritParams downloader  
#' @family downloaders
#' @seealso \url{https://github.com/axel-download-accelerator/axel/}
#' @keywords internal
axel <- function(input_url,
                 output_path,
                 background=FALSE,
                 nThread=1,
                 force_overwrite=FALSE,
                 quiet=TRUE,
                 alternate=TRUE,
                 # conda_env=NULL,
                 check_certificates=FALSE){
    message("Downloading with axel.")
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(output_path,basename(input_url))
    if(force_overwrite){
        message("+ Overwriting pre-existing file.")
        suppressWarnings(file.remove(out_file))
    }
    # axel <- CONDA.find_package("axel", conda_env = conda_env, verbose = quiet)
    axel <- "axel"
    cmd <- paste(axel,
                 input_url,
                 "-n",nThread,
                 ## Checking certificates can sometimes cause issues
                 if(check_certificates) "" else "--insecure",
                 if(force_overwrite) "" else "--no-clobber",
                 "-o",out_file,
                 if(quiet) "-q" else "",
                 # ifelse(alternate,"-a",""),
                 if(background) "& bg" else ""
    )
    # print(cmd)
    system(cmd)
    return(out_file)
}