
#' Download VCF file and its index file from Open GWAS
#' 
#' @inheritParams downloader
#' @inheritParams import_sumstats
#' @keywords internal
download_vcf <- function(vcf_url,
                         vcf_dir,
                         vcf_download=TRUE,
                         download_method="download.file",
                         force_new=FALSE,
                         quiet=TRUE,
                         nThread=1){
    #### Create save_path ####
    save_path <- file.path(vcf_dir,basename(vcf_url)) 
    
    if(file.exists(save_path) & force_new==FALSE){
        message("Using previously downloaded VCF.") 
        vcf_url <- save_path
    } else { 
        if(vcf_download){
            message("Downloading VCF ==> ",save_path)   
            dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
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
        } 
    }
    return(list(save_path=save_path,
                index_path=index_path))
}
