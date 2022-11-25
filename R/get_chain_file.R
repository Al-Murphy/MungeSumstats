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
    local_path_gunzip <- gsub(".gz", "", local_path)
    ### Download chain file ####
    if(file.exists(local_path_gunzip)){
        if(verbose)
            messager(sprintf("Using existing chain file from %s.", chain_source),v=verbose)
    }
    else{
        source_readable = ifelse(chain_source == "ensembl", "Ensembl", "UCSC Genome Browser")
        if(verbose)
            messager(sprintf("Downloading chain file from %s.", source_readable),
                     v=verbose)
        error_dwnld <-
            tryCatch(utils::download.file(remote_path, local_path),
            error = function(e){message(e);e},
            warning = function(w){message(w);w}
            )
        #if download failed use file in package
        if(is(error_dwnld,"warning")||is(error_dwnld,"error")||
                ("message" %in% names(error_dwnld) &&
                (grepl("Couldn't connect to server",error_dwnld$message)||
                    grepl("Couldn't resolve host name",error_dwnld$message)))
            ){
            chain_file <- basename(.get_ucsc_chain_remote(from, to))
            #download.file will create an empty file even if download fails
            if(file.exists(local_path)){
                rmvd <- file.remove(local_path)
                copied <- R.utils::copyFile(
                    srcPathname = system.file("extdata",chain_file,
                                              package="MungeSumstats"),
                    save_dir)
                msg <- paste0(
                    sprintf("Download of chain file from %s ", source_readable),
                    "failed, using UCSC package snapshot from 2021-10-07 ",
                    "instead.")
                message(msg)
            } 
        }
        messager(local_path,v=verbose)
        local_path <- R.utils::gunzip(local_path, overwrite=TRUE)
    }
    #### Import ####
    # Ensembl format is slightly different to UCSC, rtracklayer can't handle 
    # the spaces rather than tabs. Solution here as per RoelKluin
    # https://github.com/lawremi/rtracklayer/issues/23
    if(chain_source == "ensembl"){
        new_path = gsub(".chain", "_tabs.chain", local_path_gunzip, fixed=TRUE)
        if(!file.exists(new_path)){
            system(sprintf(
                "sed -r 's/^([0-9]+) ([0-9]+) ([0-9]+)$/\\1\\t\\2\\t\\3/' %s > %s",
                local_path_gunzip, new_path)
            )
        }
        local_path_gunzip = new_path
    }
    chain <- rtracklayer::import.chain(local_path_gunzip)
    return(chain)
}

.get_ensembl_chain_remote <- function(from = c("hg38", "hg19"),
                                      to = c("hg19", "hg38")) {
    base <- "ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/"
    ens_convert <- c("hg38" = "GRCh38", "hg19" = "GRCh37")
    ens_from <- ens_convert[from]
    ens_to <- ens_convert[to]
    return(
        paste0(base, ens_from, "_to_", ens_to, ".chain.gz")
    )
}

.get_ucsc_chain_remote <- function(from = c("hg38", "hg19"),
                                      to = c("hg19", "hg38")) {
    base <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/"
    to_caps <- paste0(toupper(substr(to, 1, 1)),substr(to, 2, nchar(to)))
    return(
        paste0(base, from, "/liftOver/", from, "To", to_caps, ".over.chain.gz")
    )
}
