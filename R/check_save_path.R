#' Check if save path is appropriate
#'
#' Returns corrected \code{save_path}, the file type, and the separator.
#' @inheritParams format_sumstats
#' @keywords internal 
check_save_path <- function(save_path,
                            write_vcf=FALSE){ 
    #### Add warning to users tha temp files arent actually saved ####
    if(dirname(save_path)==tempdir()){
        message("\n\n*****::WARNING::*****\n",
                "- Formatted results will be saved to `tempdir()` by default.\n",
                "- This means all formatted summary stats will be deleted upon ending the R session.\n",
                "- To keep formatted summary stats, change `save_path` ( e.g. `save_path=file.path('./formatted',basename(path))` ),\n",
                "  or make sure to copy files elsewhere after processing ( e.g. `file.copy(save_path, './formatted/' )`.\n",
                "*****::******::*****",
                "\n") 
    }
    
    #### Do a bit of QC to get the full path ####
    ## Expand "~" into full path bc it isn't recognized in certain envs (eg Python)
    save_path <- path.expand(save_path) 
    ## Expand relative path "./" into absolute path bc it's less ambiguous 
    save_path <- gsub("^[.]/",paste0(getwd(),"/"),save_path)
    
    suffixes <- supported_suffixes()
    if(is.null(save_path)){ 
        save_path <- paste0(tempfile(),".tsv.gz")
        file_type <- "tempfile" 
        sep <- "\t"
    } else {  
        suffix_match <- sapply(suffixes, function(x){grepl(x, tolower(save_path), ignore.case = TRUE)}) 
        if(sum(suffix_match)>0){
            file_type <- names(suffix_match)[suffix_match][1]
            sep <- if(grepl(".csv",file_type)) "," else "\t"
        } else {
            if(write_vcf==FALSE) { 
                stop("save_path file format not recognized: ",save_path,
                     "\nMust be one of: \n   ", paste(suffixes, collapse = "\n   "))
            } 
        }
    } 
    if(write_vcf){ 
        save_path <- gsub(paste(suffixes,collapse = "|"),".vcf.gz", save_path)
        sep = "\t"
        file_type = "vcf"
    } else {
        #### Account for mismatch between write_vcf and save_path ####
        suffixes.vcf <- supported_suffixes(tabular = FALSE, tabular_compressed = FALSE)
        if(any(endsWith(save_path,suffixes.vcf))){
            message("`write_vcf=FALSE` but `save_path` suggests VCF output format. ",
                    "Switching output to tabular format (.tsv.gz).")
            save_path <- gsub(paste(suffixes.vcf,collapse = "|"),".tsv.gz", save_path)
            sep = "\t"
            file_type = ".tsv"
        } 
    }  
    #### Make sure dir exists 
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE) 
    
    message("Formatted summary statistics will be saved to ==> ",save_path) 
    return(list(
        save_path=save_path,
        file_type=file_type,
        sep=sep
    ))
}
