#' Load the SNP location data for any dbSNP build
#'
#' @param ref_genome character, "GRCh37" or "GRCh38"
#' @param dbSNP integer, dbSNP build number (e.g. 144, 155, 156, …)
#' @param dbSNP_tarball optional path to a SNPlocs.Hsapiens.dbSNP<build>.GRCh<37|38>_*.tar.gz
#' @param msg optional startup message
#' @return An OnDiskLongTable of SNP locations
#' @importFrom utils untar
#' @importFrom BSgenome OnDiskLongTable
#' @export
load_snp_loc_data <- function(ref_genome,
                              dbSNP = c(144,155),
                              dbSNP_tarball = NULL,
                              msg = NULL) {
  
  # If the user has supplied a tarball, unpack it and load that namespace
  if (!is.null(dbSNP_tarball)) {
    tmpdir <- tempfile("SNPlocs_pkg")
    dir.create(tmpdir)
    utils::untar(dbSNP_tarball, exdir = tmpdir)
    # find the extracted package directory
    pkgs <- list.dirs(tmpdir, recursive = FALSE, full.names = TRUE)
    pkg_path <- pkgs[grepl(paste0("SNPlocs\\.Hsapiens\\.dbSNP", dbSNP,
                                  "\\.GRCh", ref_genome), basename(pkgs))]
    if (length(pkg_path)!=1) {
      stop("Couldn't locate SNPlocs package inside tarball")
    }
    # load its namespace  
    ns <- loadNamespace(pkg_path)
    pkg_name <- basename(pkg_path)
    SNP_LOC_DATA <- getExportedValue(ns, pkg_name)
    
  } else {
    # else fall back to the installed Bioc package
    pkg_name <- sprintf("SNPlocs.Hsapiens.dbSNP%d.GRCh%s",
                        as.integer(dbSNP), sub("GRCh","", ref_genome))
    if (!requireNamespace(pkg_name, quietly=TRUE)) {
      stop("To use dbSNP=", dbSNP,
           " on ", ref_genome,
           ", please install the Bioconductor package '",
           pkg_name, "'.")
    }
    SNP_LOC_DATA <- getExportedValue(pkg_name, pkg_name)
  }
  
  return(SNP_LOC_DATA)
}
