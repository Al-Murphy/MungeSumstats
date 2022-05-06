#' downloader wrapper
#'
#' R wrapper for 
#' \href{https://github.com/axel-download-accelerator/axel}{axel} 
#' (multi-threaded) and
#'  \link[utils]{download.file} (single-threaded) 
#'  download functions.
#' @source \href{https://stackoverflow.com/a/66602496/13214824}{
#' Suggestion to avoid 'proc$get_built_file() : Build process failed'}
#'
#' @return Local path to downloaded file.
#'
#' @param input_url input_url.
#' @param output_path output_path.
#' @param download_method \code{"axel"} (multi-threaded) or
#' \code{"download.file"} (single-threaded) .
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
                       download_method = "axel",
                       background = FALSE,
                       force_overwrite = FALSE,
                       quiet = TRUE,
                       show_progress = TRUE,
                       continue = TRUE,
                       nThread = 1,
                       alternate = TRUE,
                       check_certificates = TRUE,
                       # conda_env=NULL,
                       timeout = 30 * 60) {
    if (download_method == "axel") {
        axel_avail <- length(system("which axel", intern = TRUE)) != 0
        if (axel_avail) {
            out_file <- axel(
                input_url = input_url,
                output_path = output_path,
                background = background,
                nThread = nThread,
                force_overwrite = force_overwrite,
                quiet = quiet, # output is hella long otherwise...
                alternate = alternate,
                # conda_env=conda_env,
                check_certificates = check_certificates
            )
            if (!file.exists(out_file)) {
                message("axel download failed. Trying with download.file.")
                download_method <- "download.file"
            }
        } else {
            messager(
                "Please install axel first.\n",
                "- MacOS: brew install axel",
                "- Linux: sudo apt-get update && sudo apt-get install axel"
            )
            messager("Defaulting to download.file")
            download_method <- "download.file"
        }
    }
    if (download_method == "download.file") {
        messager("Downloading with download.file.")
        options(timeout = timeout)
        out_file <- file.path(output_path, basename(input_url))
        catch_fail <- tryCatch({
                utils::download.file(url = input_url, 
                                     destfile = out_file)
        },
            error = function(e) {message(e);e}
        )
        msg <- paste(
            "Failed to download:", input_url,
            "\nThis is likely due to either an incorrect ID/URL,",
            "or an issue with your internet connection."
        )
        #### 2nd attempt ###
        if (methods::is(catch_fail, "error") | 
            methods::is(catch_fail, "warning")) { 
            messager("Trying download.file again with different parameters.")
            catch_fail2 <- tryCatch({
                utils::download.file(url = input_url, 
                                     destfile = out_file,
                                     mode = 'wb',
                                     method = 'curl',
                                     extra = '-k')
            },
            error = function(e) {message(e);e}
            )
            
            if (methods::is(catch_fail2, "error") | 
                methods::is(catch_fail2, "warning")) {
                stop(msg)
            }
        }
        
    }
    return(out_file)
}
