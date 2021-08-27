
#' Download VCF file and its index file from Open GWAS
#'
#'
#' Ideally, we would use \href{https://github.com/MRCIEU/gwasvcf}{gwasvcf}
#' instead but it hasn't been made available on CRAN or Bioconductor yet,
#' so we can't include it as a dep.
#'
#' @return List containing the paths to the downloaded VCF and its index file.
#'
#' @param vcf_url Remote URL to VCF file.
#' @param force_new Overwrite a previously downloaded VCF with the same path name.
#' @inheritParams downloader
#' @inheritParams import_sumstats
#' @export
#' @examples
#' vcf_url <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz"
#' out_paths <- download_vcf(vcf_url = vcf_url)
download_vcf <- function(vcf_url,
                         vcf_dir = tempdir(),
                         vcf_download = TRUE,
                         download_method = "download.file",
                         force_new = FALSE,
                         quiet = TRUE,
                         nThread = 1) {
    #### Create save_path ####
    save_path <- file.path(vcf_dir, basename(vcf_url))
    index_path <- NULL

    if (file.exists(save_path) && force_new == FALSE) {
        message("Using previously downloaded VCF.")
        vcf_url <- save_path
    } else {
        if (vcf_download) {
            message("Downloading VCF ==> ", save_path)
            dir.create(dirname(save_path),
                showWarnings = FALSE,
                recursive = TRUE
            )
            #### Download main VCF file
            save_path <- downloader(
                input_url = vcf_url,
                output_path = vcf_dir,
                download_method = download_method,
                force_overwrite = force_new,
                quiet = quiet,
                nThread = nThread
            )
            #### Download tabix index file
            index_url <- paste0(vcf_url, ".tbi")
            message("Downloading VCF index ==> ", index_url)
            index_path <- downloader(
                input_url = index_url,
                output_path = vcf_dir,
                download_method = download_method,
                force_overwrite = force_new,
                quiet = quiet,
                nThread = nThread
            )
        }
    }
    return(list(
        save_path = save_path,
        index_path = index_path
    ))
}
