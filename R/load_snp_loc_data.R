#' Loads the SNP locations and alleles for Homo sapiens from dbSNP builds
#'
#' @param ref_genome   character, "GRCh37" or "GRCh38"
#' @param dbSNP        integer, dbSNP build number (144, 155, or any installed 
#' SNPlocs package)
#' @param dbSNP_tarball Optional path to a .tar.gz containing:
#'   one or more .rds files (Bioc SNPlocs package layout).
#' @param msg          optional character to message before loading
#' @return A data.table or OnDiskLongTable of SNP locations
#' @importFrom BSgenome OnDiskLongTable
#' @export
load_snp_loc_data <- function(ref_genome,
                              dbSNP,
                              dbSNP_tarball = NULL,
                              msg = NULL) {
  ## 1) Validate inputs
  ref <- toupper(ref_genome)
  if (!ref %in% c("GRCH37", "GRCH38")) {
    stop("`ref_genome` must be 'GRCh37' or 'GRCh38'.", call. = FALSE)
  }
  if (!is.null(msg)) {
    print_msg <- paste0(
      "There is no ", msg, " column found within the data. ",
      "It must be inferred from other column information."
    )
    message(print_msg)
  }
  
  ## 2) Custom‐tarball branch
  if (!is.null(dbSNP_tarball)) {
    #check ref genome matches file, else throw error
    if(!grepl(ref_genome, basename(dbSNP_tarball),ignore.case = TRUE)){
      stp_msg <- paste0(
        "The dbSNP_tarball file name does not appear to be from the specified",
        "genome build (`ref_genome`). This is based on the tarball file name.",
        "Check this is correct and update either the specified ref_genome or ",
        "the tarball file name")
      stop(stp_msg)
    }
    message("Loading SNPlocs data from tarball ", basename(dbSNP_tarball), 
            " on ", ref, ".")
    tmpdir <- tempfile("snp_loc_")
    #check if one already created - this is called multiple times per run
    tmpdirs <- list.dirs(dirname(tmpdir),recursive = FALSE)
    tmpdirs <- tmpdirs[grepl("\\bsnp_loc_",tmpdirs)]
    if(length(tmpdirs)==0){
      dir.create(tmpdir)
      utils::untar(dbSNP_tarball, exdir = tmpdir)
    }else{ #found already
      tmpdir <- tmpdirs[1] #can be multiple
    }
    # locate inst/extdata under the top‐level unpacked directory
    topdirs <- list.dirs(tmpdir, recursive = FALSE, full.names = TRUE)
    extdata <- NULL
    for (d in topdirs) {
      p <- file.path(d, "inst", "extdata")
      if (dir.exists(p)) {
        extdata <- p
      }
    }
    # fallback if unpacked flat
    if (is.null(extdata)) {
      extdata <- file.path(tmpdir, "inst", "extdata")
      if (!dir.exists(extdata)) {
        stop("Cannot locate inst/extdata in tarball", call. = FALSE)
      }
    }
    # construct and return the OnDiskLongTable
    return(BSgenome::OnDiskLongTable(extdata))
  }
  
  ## 3) Use installed SNPlocs pkg (144, 155, 156, 157, …)
  message("Loading SNPlocs data for build ", dbSNP, " on ", ref, ".")
  pkg_name <- sprintf("SNPlocs.Hsapiens.dbSNP%d.GRCh%s",
                      as.integer(dbSNP),
                      ifelse(ref == "GRCH37", "37", "38"))
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("Please install Bioconductor package '", pkg_name, "'.", call. = FALSE)
  }
  snp_loc_data <- getExportedValue(pkg_name, pkg_name)
  return(snp_loc_data)
}
