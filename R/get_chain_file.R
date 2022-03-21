#' Download chain file for liftover
#'
#' @source \href{https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/}{
#' UCSC chain files}
#' @param build_conversion converting from what build to what? hg38ToHg19 or
#' hg19ToHg38 
#' @param save_dir where is the chain file saved? Default is a temp directory
#' @param verbose extra messages printed? Default is TRUE
#' @return loaded chain file for liftover
#' @keywords internal
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @importFrom rtracklayer import.chain
get_chain_file <- function(build_conversion = c("hg38ToHg19", "hg19ToHg38"), 
                           save_dir = tempdir(),
                           verbose = TRUE) {

    #### Define paths ####
    ucsc_ref <- tolower(strsplit(build_conversion, "To")[[1]][1])
    remote_path <- file.path(
        paste0(
            "ftp://hgdownload.cse.ucsc.edu/goldenPath/", ucsc_ref,
            "/liftOver"
        ),
        paste0(build_conversion[1], ".over.chain.gz")
    )
    local_path <- file.path(save_dir, basename(remote_path))
    local_path_gunzip <- gsub(".gz", "", local_path)
    ### Download chain file ####
    if(file.exists(local_path_gunzip)){
        if(verbose)
            message("Using existing chain file.")
    }
    else{
        if(verbose)
            message("Downloading chain file from UCSC Genome Browser.")
        error_dwnld <-
            tryCatch(utils::download.file(remote_path, local_path),
            error = function(e) e,
            warning = function(w) w
            )
        #if download failed use file in package
        if(is(error_dwnld,"warning")||is(error_dwnld,"error")||
                ("message" %in% names(error_dwnld) &&
                (grepl("Couldn't connect to server",error_dwnld$message)||
                    grepl("Couldn't resolve host name",error_dwnld$message)))
            ){
            chain_file <- paste0(build_conversion[1],".over.chain.gz")
            #download.file will create an empty file even if download fails
            if(file.exists(local_path))
                rmvd <- file.remove(local_path)
            copied <- R.utils::copyFile(srcPathname =
                                   system.file("extdata",chain_file,
                                                package="MungeSumstats"),
                                   save_dir)
            msg <- paste0("Download of chain file from UCSC Genome Browser ",
                            "failed, using package snapshot from 2021-10-07 ",
                            "instead.")
            message(msg)
        }
        message(local_path)
        R.utils::gunzip(local_path, overwrite=TRUE)
    }
    #### Import ####
    chain <- rtracklayer::import.chain(local_path_gunzip)
    return(chain)
}
