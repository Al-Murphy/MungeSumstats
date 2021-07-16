

#' Import full genome-wide GWAS summary statistics from Open GWAS  
#' 
#' @param ids List of Open GWAS study IDs (e.g. \code{c("prot-a-664", "ieu-b-4760")}). 
#' @param vcf_download Download the original VCF from Open GWAS. 
#' @param vcf_dir Where to download the original VCF from Open GWAS.
#' \emph{WARNING:} This is set to \code{tempdir()} by default. 
#' This means the raw (pre-formatted) VCFs be deleted upon ending the R session.
#' Change this to keep the raw VCF file on disk (e.g. \code{vcf_dir="./raw_vcf"}). 
#' @param save_dir Directory to save formatted summary statistics in.
#' @param force_new Overwrite a previously downloaded VCF with the same path name.
#' @param parallel_across_ids If \code{parallel_across_ids=TRUE} and \code{nThread>1}, 
#' then each ID in \code{ids} will be processed in parallel. 
#' @param ... Additional arguments to be passed to \code{MungeSumstats::format_sumstats}.
#' @inheritParams format_sumstats
#' @inheritParams downloader 
#' 
#' @examples 
#' ### Search by criteria
#' metagwas <- MungeSumstats::find_sumstats(traits = c("parkinson","alzheimer"),
#'                                          min_sample_size = 5000)
#' ### Only use a subset for testing purposes                                           
#' ids <- (dplyr::arrange(metagwas, nsnp))$id[1:2]    
#' 
#' ### Default usage   
#' # datasets <- MungeSumstats::import_sumstats(ids = ids)
#'                                 
#' #### Speed up with multi-threaded download via axel
#' # datasets <- MungeSumstats::import_sumstats(ids = ids,
#' #                                            download_method="axel",
#' #                                            nThread=10,
#' #                                            parallel_across_ids=TRUE)
#' @return Either a named list of data objects or paths, 
#' depending on the arguments passed to \code{format_sumstats}.
#' @export 
import_sumstats <- function(ids,  
                            vcf_dir=tempdir(),
                            vcf_download=TRUE,
                            save_dir=tempdir(),
                            write_vcf=FALSE,
                            download_method="download.file",
                            quiet=TRUE, 
                            force_new=FALSE,
                            nThread=1, 
                            parallel_across_ids=FALSE,
                            ...){  
    # vcf_dir=tempdir(); vcf_download=TRUE;download_method="axel";quiet=FALSE;force_new=FALSE;nThread=10; ids=c("ieu-a-1124","ieu-a-1125"); id=ids[1]; ref_genome=NULL;
    start_all <- Sys.time()
    ids <- unique(ids)
    message("Processing ",length(ids)," datasets from Open GWAS.")   
    #### Handle parallelization parameters ####
    if(parallel_across_ids) {
        nThread_acrossIDs <- nThread
        nThread <- 1
        message("`parallel_across_ids=TRUE` : most messages will be hidden in this mode.")
    } else {nThread_acrossIDs <- 1}
    
    ouputs <- parallel::mclapply(ids, function(id){ 
        out <-  tryCatch(expr = {
            start <- Sys.time()
            message_parallel("\n========== Processing dataset : ",id," ==========\n") 
            vcf_url <- file.path("https://gwas.mrcieu.ac.uk/files",id,paste0(id,".vcf.gz")) 
            #### Optional:: download VCF ####
            vcf_paths <- download_vcf(vcf_url=vcf_url,
                                      vcf_dir=vcf_dir, 
                                      vcf_download=vcf_download,
                                      download_method=download_method,
                                      force_new=force_new,
                                      quiet=quiet,
                                      nThread=nThread) 
            #### format_sumstats ####
            save_path <- file.path(save_dir,basename(vcf_url))
            if(!write_vcf) save_path <- gsub(".vcf.gz",".tsv.gz",save_path)
            reformatted <- format_sumstats(path=vcf_paths$save_path,
                                           save_path=save_path,
                                           force_new=force_new,
                                           nThread=nThread, 
                                           ...)   
            end <- Sys.time() 
            message("\n",id," : Done in ",round(difftime(end, start, units='mins'), 2)," minutes.")
            return(reformatted) 
        }, error = function(e){message(e);return(as.character(e))}) 
         
        return(out)
    }, mc.cores = nThread_acrossIDs) %>% `names<-`(ids)
    
    end_all <- Sys.time() 
    message("\nDone with all processing in ",
            round(difftime(end_all, start_all, units='mins'), 2)," minutes.")
    return(ouputs) 
}