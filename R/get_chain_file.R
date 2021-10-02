#' Download chain file for liftover
#'
#' @source \href{https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/}{
#' UCSC chain files}
#' @param build_conversion converting from what build to what? hg38ToHg19 or
#' hg19ToHg38
#' @param ucsc_ref converting from what? hg19 or hg38
#' @param save_dir where is the chain file saved? Default is a temp directory
#' @param verbose extra messages printed? Default is TRUE
#' @return loaded chain file for liftover
#' @keywords internal
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @importFrom rtracklayer import.chain
get_chain_file <- function(build_conversion = c("hg38ToHg19", "hg19ToHg38"),
                           ucsc_ref = c("hg19", "hg38"),
                           save_dir = tempdir(),
                           verbose = TRUE) {

    #### Define paths ####
    remote_path <- file.path(
        paste0(
            "ftp://hgdownload.cse.ucsc.edu/goldenPath/", ucsc_ref[1],
            "/liftOver"
        ),
        paste0(build_conversion[1], ".over.chain.gz")
    )
    local_path <- file.path(save_dir, basename(remote_path))
    local_path_gunzip <- gsub(".gz", "", local_path)
    ### Download chain file ####
    if (file.exists(local_path_gunzip)) {
        if (verbose) {
              message("Using existing chain file.")
          }
    } else {
        if (verbose) {
              message("Downloading chain file from UCSC Genome Browser.")
          }
        utils::download.file(remote_path, local_path)
        print(local_path)
        R.utils::gunzip(local_path)
    }
    #### Import ####
    chain <- rtracklayer::import.chain(local_path_gunzip)
    return(chain)
}
