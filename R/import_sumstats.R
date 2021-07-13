

#' Import full genome-wide GWAS summary statistics from Open GWAS  
#' 
#' @param ids List of Open GWAS study IDs (e.g. \code{c("prot-a-664", "ieu-b-4760")}). 
#' @param vcf_download Download the original VCF from Open GWAS. 
#' @param vcf_dir Where to download the original VCF from Open GWAS.
#' @param save_dir Directory to save formatted summary statistics in.
#' @param force_new Overwrite a previously downloaded VCF with the same path name.
#' @param ... Additional arguments to be passed to \code{MungeSumstats::format_sumstats}.
#' @inheritParams format_sumstats
#' @inheritParams downloader 
#' 
#' @examples 
#' ### Search by criteria
#' metagwas <- MungeSumstats::find_sumstats()
#' ### Only use a subset for testing purposes                                           
#' ids <- (dplyr::arrange(metagwas, nsnp))$id[1:2]    
#' 
#' ### Default usage   
#' # datasets <- MungeSumstats::import_sumstats(ids = ids)
#'                                 
#' #### Speed up with multi-threaded download via axel 
#' # datasets <- MungeSumstats::import_sumstats(ids = ids,  download_method="axel", nThread=10)                                           
#' 
#' @return Either a named list of data objects or paths, depending on the arguments passed to 
#' @export
#' @importFrom VariantAnnotation readVcf geno 
import_sumstats <- function(ids,  
                            vcf_dir=tempdir(),
                            vcf_download=TRUE,
                            save_dir="./formatted",
                            write_vcf=FALSE,
                            download_method="download.file",
                            quiet=TRUE, 
                            force_new=FALSE,
                            nThread=1, 
                            ...){  
    # vcf_dir=tempdir(); vcf_download=TRUE;download_method="axel";quiet=FALSE;force_new=FALSE;nThread=10; ids=c("ieu-a-1124","ieu-a-1125"); id=ids[1]; ref_genome=NULL;
    
    ids <- unique(ids)
    message("Processing ",length(ids)," datasets from Open GWAS.") 
    
    ouputs <- lapply(ids, function(id){
        start <- Sys.time()
        out <-  tryCatch(expr = {
            message("\n========== Processing dataset : ",id," ==========\n") 
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
            return(reformatted) 
        }, error = function(e){message(e);return(as.character(e))}) 
        
        end <- Sys.time() 
        message("\n",id," : Done in ",round(difftime(end, start, units='mins'), 2)," minutes.")
        return(out)
    }) %>% `names<-`(ids)
    return(ouputs) 
}