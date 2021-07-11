

#' Import full genome-wide GWAS summary statistics from Open GWAS  
#' 
#' @param ids List of Open GWAS study IDs (e.g. \code{c("prot-a-664", "ieu-b-4760")}). 
#' @param vcf_download Download the original VCF from Open GWAS. 
#' @param vcf_dir Where to download the original VCF from Open GWAS.
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
#' datasets <- MungeSumstats::import_sumstats(ids = ids)                                           
#' @return Either a named list of data objects or paths, depending on the arguments passed to 
#' @export
#' @importFrom VariantAnnotation readVcf geno 
import_sumstats <- function(ids, 
                            vcf_download=FALSE,
                            vcf_dir=tempdir(),
                            download_method="download.file",
                            quiet=TRUE,
                            force_new=FALSE,
                            nThread=1, 
                            ...){  
    ids <- unique(ids)
    message("Processing ",length(ids)," datasets from Open GWAS.") 
    
    ouputs <- lapply(ids, function(id){
        start <- Sys.time()
        out <-  tryCatch(expr = {
            message("\n========== Processing dataset : ",id," ==========\n") 
            vcf_url <- file.path("https://gwas.mrcieu.ac.uk/files",id,paste0(id,".vcf.gz"))
            
            #### Create save_path ####
            save_path <- file.path(vcf_dir,basename(vcf_url)) 
            dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
            #### Optional:: download VCF ####
            if(file.exists(save_path) & force_new==FALSE){
                message("Using previously downloaded VCF.") 
                vcf_url <- save_path
            } else { 
                if(vcf_download){
                    message("Downloading VCF ==> ",save_path)  
                    #### Download main VCF file
                    save_path <- downloader(input_url = vcf_url, 
                                            output_path = vcf_dir, 
                                            download_method = download_method,
                                            force_overwrite = force_new,
                                            quiet = quiet,
                                            nThread = nThread)
                    #### Download tabix index file
                    index_path <- downloader(input_url = paste0(vcf_url,".tbi"), 
                                             output_path = paste0(save_path,".tbi"), 
                                             download_method = download_method,
                                             force_overwrite = force_new,
                                             quiet = quiet,
                                             nThread = nThread)
                    vcf_url <- save_path 
                } 
            }
            reformatted <- format_sumstats(path=vcf_url,
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