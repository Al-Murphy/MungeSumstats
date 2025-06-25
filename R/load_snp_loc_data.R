#' @title Load SNP location data for any dbSNP build
#' @description
#' Loads the SNP locations and alleles for Homo sapiens from dbSNP builds.  
#' For builds distributed as Bioconductor packages (e.g. 144, 155) it will  
#' load the installed SNPlocs.Hsapiens.dbSNP<build>.GRCh<37|38> object.  
#' For newer builds provided as Bioc‐style source bundles (tar.gz of a  
#' SNPlocs package), pass the path via `dbSNP_tarball` to load directly.
#'
#' @param ref_genome    character, "GRCh37" or "GRCh38"
#' @param dbSNP         integer, dbSNP build number (e.g. 144,155,156,157,…)
#' @param dbSNP_tarball optional path to a SNPlocs.Hsapiens.dbSNP<build>.GRCh<37|38>_*.tar.gz
#' @param msg           optional startup message
#' @return An OnDiskLongTable of SNP locations
#' @importFrom utils untar
#' @importFrom BSgenome OnDiskLongTable
#' @export
load_snp_loc_data <- function(ref_genome,
                              dbSNP,
                              dbSNP_tarball = NULL,
                              msg = NULL) {
  ## 1) Validate inputs
  ref <- toupper(ref_genome)
  if (!ref %in% c("GRCH37","GRCH38")) {
    stop("`ref_genome` must be 'GRCh37' or 'GRCh38'.", call. = FALSE)
  }
  if (!is.null(msg)) message(msg)
  message("Loading SNPlocs data for dbSNP build ", dbSNP,
          " on reference ", ref, ".")
  
  ## 2) Tarball branch: user‐supplied Bioc‐style SNPlocs source bundle
  if (!is.null(dbSNP_tarball)) {
    if (!file.exists(dbSNP_tarball)) {
      stop("Tarball not found: ", dbSNP_tarball, call. = FALSE)
    }
    tmpdir <- tempfile("snp_loc_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    utils::untar(dbSNP_tarball, exdir = tmpdir)
    
    # locate inst/extdata under the top‐level unpacked directory
    topdirs <- list.dirs(tmpdir, recursive = FALSE, full.names = TRUE)
    extdata <- NULL
    for (d in topdirs) {
      p <- file.path(d, "inst", "extdata")
      if (dir.exists(p)) {
        extdata <- p
        break
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
  
  ## 3) Installed‐package branch: SNPlocs.Hsapiens.dbSNP<build>.GRCh<37|38>
  pkg_name <- sprintf("SNPlocs.Hsapiens.dbSNP%d.GRCh%s",
                      as.integer(dbSNP),
                      if (ref == "GRCH37") "37" else "38")
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("No tarball supplied and cannot load package '", pkg_name,
         "'. Please install it or provide `dbSNP_tarball`.", call. = FALSE)
  }
  
  required <- c("breakpoints.rds","spatial_index.rds","header.rds","seqinfo.txt")
  found    <- basename(list.files(extdata))
  missing  <- setdiff(required, found)
  if (length(missing)) {
    stop("Tarball is missing: ", paste(missing, collapse=", "), call.=FALSE)
  }
  # get the OnDiskLongTable object exported by the SNPlocs package
  return(utils:::getExportedValue(pkg_name, pkg_name))
}
