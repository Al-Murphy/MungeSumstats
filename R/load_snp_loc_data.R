#' Loads the SNP locations and alleles for Homo sapiens from dbSNP builds
#'
#' @param ref_genome   character, "GRCh37" or "GRCh38"
#' @param dbSNP        integer, dbSNP build number (144, 155, or any installed SNPlocs package)
#' @param dbSNP_tarball Optional path to a .tar.gz containing:
#'   one or more .rds files (Bioc SNPlocs package layout).
#' @param msg          optional character to message before loading
#' @return A data.table or OnDiskLongTable of SNP locations
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
    message(msg)
  }
  message("Loading SNPlocs data for build ", dbSNP, " on ", ref, ".")
  
  ## 2) Custom‐tarball branch
  if (!is.null(dbSNP_tarball)) {
    tmpdir <- tempfile("snp_loc_"); dir.create(tmpdir)
    utils::untar(dbSNP_tarball, exdir = tmpdir)
    
    ## 2a) Flat TSV?
    tsvs <- list.files(tmpdir, "\\.tsv$", full.names = TRUE, recursive = TRUE)
    if (length(tsvs)) {
      return(data.table::fread(tsvs[[1]]))
    }
    
    ## 2b) RDS blocks?
    rds_files <- list.files(tmpdir, "\\.rds$", full.names = TRUE, recursive = TRUE)
    if (length(rds_files)) {
      dt_list <- lapply(rds_files, readRDS)
      return(data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE))
    }
    
    stop("No .tsv or .rds found in tarball: ", dbSNP_tarball, call. = FALSE)
  }
  
  ## 3) Fallback to installed SNPlocs pkg (144, 155, 156, 157, …)
  pkg_name <- sprintf("SNPlocs.Hsapiens.dbSNP%d.GRCh%s",
                      as.integer(dbSNP),
                      ifelse(ref == "GRCH37", "37", "38"))
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("Please install Bioconductor package '", pkg_name, "'.", call. = FALSE)
  }
  snp_loc_data <- getExportedValue(pkg_name, pkg_name)
  return(snp_loc_data)
}
