

#' downloader wrapper
#'
#' R wrapper for \code{"axel"} (multi-threaded) and
#'  \code{"download.file"} (single-threaded) download functions. 
#'  
#' @return Local path to downloaded file. 
#' 
#' @param input_url input_url.
#' @param output_path output_path.
#' @param download_method \code{"axel"} (multi-threaded) or \code{"download.file"} (single-threaded) .
#' @param background Run in background 
#' @param force_overwrite Overwrite existing file.
#' @param quiet Run quietly.
#' @param show_progress show_progress.
#' @param continue continue.
#' @param nThread Number of threads to parallelize over.
#' @param alternate alternate,
#' @param check_certificates check_certificates
#' @param timeout How many seconds before giving up on download. 
#' Passed to \code{download.file}. Default: \code{30*60} (30min). 
#' @family downloaders
#' @keywords internal
downloader <- function(input_url,
                       output_path,
                       download_method="axel",
                       background=FALSE,
                       force_overwrite=FALSE,
                       quiet=TRUE,
                       show_progress=TRUE,
                       continue=TRUE, 
                       nThread=1,
                       alternate=TRUE,
                       check_certificates=TRUE,
                       # conda_env=NULL,
                       timeout=30*60){
    if(download_method=="axel"){
        axel_avail <- length(system("which axel",intern = TRUE))!=0
        if(axel_avail 
           # | !is.null(conda_env)
           ){
            out_file <- axel(input_url=input_url,
                             output_path=output_path,
                             background=background,
                             nThread=nThread,
                             force_overwrite=force_overwrite,
                             quiet=quiet, # output is hella long otherwise...
                             alternate=alternate,
                             # conda_env=conda_env,
                             check_certificates=check_certificates)
            if(!file.exists(out_file)){
                message("axel download failed. Trying with download.file.")
                download_method <- "download.file"
            } 
        } else {
            message("axel not installed.\n",
                    "For Mac users, please install via brew in the command line (`brew install axel`)");
            message("Defaulting to download.file")
            download_method <- "download.file"
        }
        
    }
    if(download_method=="download.file"){ 
        message("Downloading with download.file.")
        options(timeout=timeout) 
        out_file <- file.path(output_path, basename(input_url))
        utils::download.file(input_url, out_file) 
    }
    return(out_file)
}
 
