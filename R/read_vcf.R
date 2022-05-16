#' Read in VCF file
#' 
#' Read in a VCF file as a \link[VariantAnnotation]{VCF} or a 
#' \link[data.table]{data.table}. 
#' Can optionally save the VCF/data.table as well. 
#' @inheritParams format_sumstats
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
#' @param verbose Print messages.
#' @inheritParams check_empty_cols
#' @inheritParams VariantAnnotation::ScanVcfParam
#'
#' @return The VCF file in data.table format.
#' @export
#' @importFrom VariantAnnotation readVcf writeVcf vcfSamples
#' @importFrom GenomicRanges seqinfo
#' @importFrom data.table as.data.table setnames fwrite 
#' @importFrom methods as show
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
#' @examples 
#' #### Local file ####
#' path <- system.file("extdata","ALSvcf.vcf", package="MungeSumstats")
#' sumstats_dt <- read_vcf(path = path)
#' 
#' #### Remote file ####
#' path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz" 
#' sumstats_dt2 <- read_vcf(path = path)
read_vcf <- function(path, 
                     write_vcf = FALSE,
                     save_path = NULL,
                     tabix_index = FALSE,
                     samples = 1,
                     which = NULL,
                     use_params = TRUE,
                     sampled_rows = 1e7,
                     nThread = 1,
                     verbose = TRUE){
    #### Read #### 
    {
        t1 <- Sys.time()
        if((!is.null(which)) || isTRUE(use_params)){
            #### Make sure file is compressed and indexed ####
            ## File must be indexed in order to use param 
            ## (even if only specifying columns) 
            path <- index_vcf(path = path,
                              verbose = verbose)
            param <- select_vcf_fields(path = path, 
                                       which = which, 
                                       samples = samples,
                                       sampled_rows = sampled_rows,
                                       nThread = nThread)
        } else {
            param <- VariantAnnotation::ScanVcfParam()
        }
        messager("Reading VCF file.")
        vcf <- suppressWarnings(
            VariantAnnotation::readVcf(file = path,
                                       param = param)
        )
        methods::show(round(difftime(Sys.time(),t1),1))
    } 
    #### Check save path ####
    if(!is.null(save_path)){
        check_save_path_out <- check_save_path(
            ### dummy arg, not used here
            log_folder = "NULL", 
            ### dummy arg, not used here
            log_folder_ind = FALSE, 
            save_path = save_path, 
            write_vcf = write_vcf,
            tabix_index = tabix_index)
        save_path <- check_save_path_out$save_path
    }
    #### Save as VCF ####
    if(isTRUE(write_vcf)) {
        if(!is.null(save_path)){ 
            VariantAnnotation::writeVcf(
                obj = vcf,
                filename = check_save_path_out$save_path,
                index = tabix_index
            )
        }
        messager("Returning as VCF.")
        return(vcf)
    }
    #### Convert to data.table ####  
    samples <- VariantAnnotation::vcfSamples(param)
    sumstats_dt <- vcf2df(vcf = vcf, 
                          add_sample_names = length(samples)!=1)
    sample_id <- rownames(vcf@colData) 
    remove(vcf) 
    #### Remove duplicated columns ####
    drop_duplicate_cols(dt = sumstats_dt)
    #### Remove sample suffix ####
    data.table::setnames(
        x = sumstats_dt,
        old = colnames(sumstats_dt), 
        new = gsub(paste0("_",sample_id,collapse = "|"),"",
                   colnames(sumstats_dt), 
                   ignore.case = TRUE)
    )   
    #### Remove empty columns #####
    sumstats_dt <- remove_empty_cols(sumstats_dt = sumstats_dt,
                                     sampled_rows = sampled_rows)
    #### Unlist columns inplace ####
    unlist_dt(dt = sumstats_dt) 
    #### Remove duplicated rows #### 
    sumstats_dt <- unique(sumstats_dt)
    #### Prepare SNP column ####
    read_vcf_markername(sumstats_dt = sumstats_dt) 
    #### Prepare/un-log P col ####
    read_vcf_pval(sumstats_dt = sumstats_dt)  
    #### Prepare INFO col ####
    read_vcf_info(sumstats_dt = sumstats_dt)  
    #### Rename start col #####
    data.table::setnames(sumstats_dt,"start","BP")
    #### Write new data ####
    if (!is.null(save_path)) { 
        messager("Storing intermediate file before proceeding ==>",save_path)
        data.table::fwrite(
            x = sumstats_dt,
            file = save_path, 
            sep = "\t",
            nThread = nThread,
        )
    }
    return(sumstats_dt)
}