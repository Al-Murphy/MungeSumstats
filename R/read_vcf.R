#' Read in VCF file
#' 
#' Read in a VCF file as a \link[VariantAnnotation]{VCF} or a 
#' \link[data.table]{data.table}. 
#' Can optionally save the VCF/data.table as well. 
#' 
#' @param path Path to local or remote VCF file.
#' @param as_datatable Return the data as a
#'  \link[data.table]{data.table} (default: \code{TRUE}) 
#'  or a \link[VariantAnnotation]{VCF} (\code{FALSE}).
#' @param samples Which samples to use:
#' \itemize{
#' \item{1 : }{Only the first sample will be used (\emph{DEFAULT}).}
#' \item{NULL : }{All samples will be used.}
#' \item{c("<sample_id1>","<sample_id2>",...) : }{
#' Only user-selected samples will be used (case-insensitive).}
#' }
#' @param use_params 
#' When \code{TRUE} (default), increases the speed of reading in the VCF by
#' omitting columns that are empty based on the head of the VCF (NAs only). 
#' NOTE that that this requires the VCF to be sorted, bgzip-compressed, 
#' tabix-indexed, which \link[MungeSumstats]{read_vcf} will attempt to do.
#' @param download Download the VCF (and its index file) 
#' to a temp folder before reading it into R. 
#' This is important to keep \code{TRUE} when \code{nThread>1} to avoid
#' making too many queries to remote file. 
#' @param mt_thresh When the number of rows (variants) in the VCF is 
#' \code{< mt_thresh}, only use single-threading for reading in the VCF.
#' This is because the overhead of parallelisation outweighs the speed benefits
#' when VCFs are small. 
#' @param verbose Print messages.
#' @inheritParams check_empty_cols
#' @inheritParams format_sumstats
#' @inheritParams downloader
#' @inheritParams download_vcf
#' @inheritParams import_sumstats
#' @inheritParams VariantAnnotation::ScanVcfParam
#'
#' @return The VCF file in data.table format.
#' @export
#' @importFrom VariantAnnotation writeVcf 
#' @importFrom data.table as.data.table setnames fwrite 
#' @source 
#' \code{
#' #### Benchmarking ####
#' library(VCFWrenchR)
#' library(VariantAnnotation)
#' path <- "https://gwas.mrcieu.ac.uk/files/ubm-a-2929/ubm-a-2929.vcf.gz"
#' vcf <- VariantAnnotation::readVcf(file = path)
#' N <- 1e5
#' vcf_sub <- vcf[1:N,]
#' res <- microbenchmark::microbenchmark(
#'     "vcf2df"={dat1 <- MungeSumstats:::vcf2df(vcf = vcf_sub)},
#'     "VCFWrenchR"= {dat2 <- as.data.frame(x = vcf_sub)},
#'     "VRanges"={dat3 <- data.table::as.data.table(
#'         methods::as(vcf_sub, "VRanges"))},
#'     times=1
#' )
#' }
#' @source \href{https://github.com/Bioconductor/VariantAnnotation/issues/57}{
#' Discussion on VariantAnnotation GitHub}
#' @source \href{https://github.com/Bioconductor/VariantAnnotation/issues/59}{
#' Discussion on VariantAnnotation GitHub} 
#' @examples 
#' #### Local file ####
#' path <- system.file("extdata","ALSvcf.vcf", package="MungeSumstats")
#' sumstats_dt <- read_vcf(path = path)
#' 
#' #### Remote file ####
#' ## Small GWAS (0.2Mb)
#' # path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz" 
#' # sumstats_dt2 <- read_vcf(path = path)
#' 
#' ## Large GWAS (250Mb)
#' # path <- "https://gwas.mrcieu.ac.uk/files/ubm-a-2929/ubm-a-2929.vcf.gz"
#' # sumstats_dt3 <- read_vcf(path = path, nThread=11)
#' 
#' ### Very large GWAS (500Mb)
#' # path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-1124/ieu-a-1124.vcf.gz"
#' # sumstats_dt4 <- read_vcf(path = path, nThread=11)
read_vcf <- function(path, 
                     as_datatable = TRUE,
                     save_path = NULL,
                     tabix_index = FALSE,
                     samples = 1,
                     which = NULL,
                     use_params = TRUE,
                     sampled_rows = 1e4L,
                     download = TRUE,
                     vcf_dir = tempdir(),
                     download_method = "download.file",
                     force_new = FALSE,
                     mt_thresh = 1e5L,
                     nThread = 1,
                     verbose = TRUE){ 
    #### Read VCF #### 
    ## Returns either a VCF or a data.table, depending on as_datatable arg.
    vcf_dt <- read_vcf_parallel(path = path,
                                samples = samples,
                                which = which,
                                use_params = use_params,
                                as_datatable = as_datatable,
                                sampled_rows = sampled_rows,
                                
                                download = download,
                                vcf_dir = vcf_dir,
                                download_method = download_method,
                                force_new = force_new, 
                                nThread = nThread,
                                verbose = verbose)
    #### Check save path ####
    if(!is.null(save_path)){
        check_save_path_out <- check_save_path(
            ### dummy arg, not used here
            log_folder = "NULL", 
            ### dummy arg, not used here
            log_folder_ind = FALSE, 
            save_path = save_path, 
            write_vcf = !as_datatable,
            tabix_index = tabix_index)
        save_path <- check_save_path_out$save_path
    }
    #### Save as VCF ####
    if(isFALSE(as_datatable)) {
        if(!is.null(save_path)){ 
            VariantAnnotation::writeVcf(
                obj = vcf_dt,
                filename = check_save_path_out$save_path,
                index = tabix_index
            )
        }
        messager("Returning as VCF.")
        return(vcf_dt)
    } 
    #### Prepare SNP column ####
    read_vcf_markername(sumstats_dt = vcf_dt) 
    #### Prepare/un-log P col ####
    read_vcf_pval(sumstats_dt = vcf_dt)  
    #### Prepare INFO col ####
    read_vcf_info(sumstats_dt = vcf_dt)  
    #### Rename start col #####
    data.table::setnames(vcf_dt,"start","BP")
    #### Write new data ####
    if (!is.null(save_path)) { 
        messager("Storing intermediate file before proceeding ==>",save_path,
                 v=verbose)
        data.table::fwrite(
            x = vcf_dt,
            file = save_path, 
            sep = "\t",
            nThread = nThread,
            na = "NA",
            quote = FALSE
        )
    } 
    return(vcf_dt)
}