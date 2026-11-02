#' Download chain file for liftover
#'
#' @source \href{https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/}{
#' UCSC chain files}
#' @source \href{https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/}{
#' Ensembl chain files}
#' @param from genome build converted from ("hg38", "hg19")
#' @param to genome build converted to ("hg19", "hg38")
#' @param chain_source chain file source used ("ucsc" as default, or "ensembl")
#' @param save_dir where is the chain file saved? Default is a temp directory
#' @param verbose extra messages printed? Default is TRUE
#' @return loaded chain file for liftover
#' @keywords internal
#' @importFrom utils download.file
#' @importFrom R.utils gunzip copyFile
#' @importFrom rtracklayer import.chain
get_chain_file <- function(from = c("hg38", "hg19"),
                           to = c("hg19", "hg38"),
                           chain_source = c("ucsc", "ensembl"),
                           save_dir = tempdir(),
                           verbose = TRUE) {
  #### Match args ####
  from <- match.arg(from)
  to <- match.arg(to)
  chain_source <- match.arg(chain_source)
  if(from == to){
    stop("Cannot get a chain file from one reference to the same.")
  }
  if(chain_source == "ucsc"){
    message("Note that you are fetching the UCSC chain file, ",
            "which requires a licence for commercial use.")
  }
  
  #### Define paths ####
  remote_path = switch(
    chain_source,
    "ensembl" = .get_ensembl_chain_remote(from, to),
    "ucsc" = .get_ucsc_chain_remote(from, to)
  )
  local_path <- file.path(save_dir, basename(remote_path))
  local_path_gunzip <- sub("\\.gz$", "", local_path)
  ### Download chain file ####
  if(file.exists(local_path_gunzip)){
    if(verbose)
      messager(sprintf("Using existing chain file from %s.", 
                       chain_source),v=verbose)
  }
  else{
    source_readable = ifelse(chain_source == "ensembl", "Ensembl", 
                             "UCSC Genome Browser")
    if(verbose)
      messager(sprintf("Downloading chain file from %s.", 
                       source_readable),
               v=verbose)
    error_dwnld <-
      tryCatch(utils::download.file(remote_path, local_path),
               error = function(e){e},
               warning = function(w){w}
      )
    #if download failed use file in package
    if(!(is.numeric(error_dwnld) && error_dwnld == 0) ||
       !file.exists(local_path)){
      #a failed download may or may not leave an empty file behind
      if(file.exists(local_path))
        rmvd <- file.remove(local_path)
      if(chain_source != "ensembl"){
        stop(sprintf("Download of chain file from %s ", source_readable),
             "failed, use ensembl chain file (chain_source='ensembl') ",
             "\nfor an offline solution.")
      }
      copied <- file.copy(
        system.file("extdata", basename(remote_path),
                    package = "MungeSumstats"),
        local_path)
      messager(paste0(
        sprintf("Download of chain file from %s ", source_readable),
        "failed, using ensembl chain file downloaded on 2022-11-25 ",
        "instead."))
    }
    messager(local_path,v=verbose)
    local_path <- R.utils::gunzip(local_path, overwrite=TRUE)
    local_path_gunzip <- local_path
  }
  #### Import ####
  #NOT NECESSARY WITH FIX TO rtracklayer - v1.59.1
  #(https://github.com/lawremi/rtracklayer/issues/23
  # Ensembl format is slightly different to UCSC, rtracklayer can't handle 
  # the spaces rather than tabs. Solution here as per RoelKluin
  # https://github.com/lawremi/rtracklayer/issues/23
  #if(chain_source == "ensembl"){
  #    new_path = gsub(".chain", "_tabs.chain", local_path_gunzip, fixed=TRUE)
  #    if(!file.exists(new_path)){
  #        system(sprintf(
  #            "sed -E 's/^([0-9]+) ([0-9]+) ([0-9]+)$/\\1\\t\\2\\t\\3/' %s > %s",
  #            local_path_gunzip, new_path)
  #        )
  #    }
  #    local_path_gunzip <- new_path
  #}
  chain <- rtracklayer::import.chain(local_path_gunzip)
  return(chain)
}

.get_ensembl_chain_remote <- function(from = c("hg38", "hg19"),
                                      to = c("hg19", "hg38")) {
  base <- "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/"
  ens_convert <- c("hg38" = "GRCh38", "hg19" = "GRCh37")
  ens_from <- ens_convert[from]
  ens_to <- ens_convert[to]
  return(
    paste0(base, ens_from, "_to_", ens_to, ".chain.gz")
  )
}

.get_ucsc_chain_remote <- function(from = c("hg38", "hg19"),
                                   to = c("hg19", "hg38")) {
  base <- "https://hgdownload.soe.ucsc.edu/goldenPath/"
  to_caps <- paste0(toupper(substr(to, 1, 1)),substr(to, 2, nchar(to)))
  return(
    paste0(base, from, "/liftOver/", from, "To", to_caps, ".over.chain.gz")
  )
}
