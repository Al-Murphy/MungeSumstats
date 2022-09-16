#' Import full genome-wide GWAS summary statistics from Open GWAS
#'
#' Requires internet access to run.
#'
#' @param ids List of Open GWAS study IDs 
#' (e.g. \code{c("prot-a-664", "ieu-b-4760")}).
#' @param vcf_download Download the original VCF from Open GWAS.
#' @param vcf_dir Where to download the original VCF from Open GWAS.
#' \emph{WARNING:} This is set to \code{tempdir()} by default.
#' This means the raw (pre-formatted) VCFs be deleted upon ending the R session.
#' Change this to keep the raw VCF file on disk 
#' (e.g. \code{vcf_dir="./raw_vcf"}).
#' @param save_dir Directory to save formatted summary statistics in.
#' @param force_new_vcf Overwrite a previously downloaded VCF 
#' with the same path name.
#' @param parallel_across_ids If \code{parallel_across_ids=TRUE} 
#' and \code{nThread>1}, 
#' then each ID in \code{ids} will be processed in parallel.
#' @param ... Additional arguments passed to 
#' \link[MungeSumstats]{format_sumstats}.
#' @inheritParams format_sumstats
#' @inheritDotParams format_sumstats
#' @inheritParams downloader
#'
#' @return Either a named list of data objects or paths,
#' depending on the arguments passed to \code{format_sumstats}.
#'
#' @examples
#' #only run the examples if user has internet access:
#' if(try(is.character(getURL("www.google.com")))==TRUE){
#' ### Search by criteria
#' metagwas <- find_sumstats(
#'     traits = c("parkinson", "alzheimer"),
#'     min_sample_size = 5000
#' )
#' ### Only use a subset for testing purposes
#' ids <- (dplyr::arrange(metagwas, nsnp))$id
#'
#' ### Default usage
#' ## You can supply \code{import_sumstats()}
#' ## with a list of as many OpenGWAS IDs as you want,
#' ## but we'll just give one to save time.
#'
#' ## Call uses reference genome as default with more than 2GB of memory,
#' ## which is more than what 32-bit Windows can handle so remove certain checks
#' ## commented out down to runtime
#' # datasets <- import_sumstats(ids = ids[1])
#' }
#' @export
#' @importFrom dplyr %>%
#' @importFrom RCurl getURL
import_sumstats <- function(ids,
                            vcf_dir = tempdir(),
                            vcf_download = TRUE,
                            save_dir = tempdir(),
                            write_vcf = FALSE,
                            download_method = "download.file",
                            quiet = TRUE,
                            force_new = FALSE,
                            force_new_vcf = FALSE,
                            nThread = 1,
                            parallel_across_ids = FALSE,
                            ...) {
    # vcf_dir=tempdir(); vcf_download=TRUE;download_method="axel";quiet=FALSE;
    # force_new=FALSE;nThread=10; ids=c("ieu-a-1124","ieu-a-1125"); id=ids[1];
    # ref_genome=NULL;
    
    #### Check download method ####
    msg_dwnld <- paste0(
        "download_method must be `download.file` (single-threaded)",
        "or `axel` (multi-threaded)."
    )
    if (!tolower(download_method) %in% c("download.file", "axel")) {
        stop(msg_dwnld)
    }
    download_method <- tolower(download_method)[1] 
    #### start overall timer ####
    start_all <- Sys.time()
    ids <- unique(ids)
    message("Processing ", length(ids), " datasets from Open GWAS.")
    #### Handle parallelization parameters ####
    if (parallel_across_ids) {
        nThread_acrossIDs <- nThread
        nThread <- 1
        message(
            "parallel_across_ids=TRUE: ",
            "most notes will be hidden in this mode."
        )
    } else {
        nThread_acrossIDs <- 1
    }
    #### Iterate over IDS ####
    ouputs <- parallel::mclapply(ids, function(id) {
        out <- tryCatch(expr = {
            start <- Sys.time()
            message_parallel(
                "\n========== Processing dataset : ", id,
                " ==========\n"
            )
            vcf_url <-
                paste(
                    "https://gwas.mrcieu.ac.uk/files", id,
                    paste0(id,".vcf.gz"),sep="/"
                )
            #### Create path of the output file ####
            id_dir <- file.path(save_dir, id)
            save_path <- file.path(id_dir, basename(vcf_url))
            if (!write_vcf) save_path <- gsub(".vcf.gz", ".tsv.gz", save_path)
            dir.create(id_dir, 
                       showWarnings = FALSE, recursive = TRUE)
            #### Create dataset-specific folders/paths ####
            log_folder <- file.path(id_dir,"logs")
            dir.create(log_folder, showWarnings = FALSE, recursive = TRUE) 
            #### Check if formatted file exists BEFORE downloading VCF ####
            if((!file.exists(save_path)) || isTRUE(force_new)){  
                #### Optional:: download VCF ####
                vcf_paths <- download_vcf(
                    vcf_url = vcf_url,
                    vcf_dir = vcf_dir,
                    vcf_download = vcf_download,
                    download_method = download_method,
                    force_new = force_new_vcf,
                    quiet = quiet,
                    nThread = nThread
                ) 
            } 
            #### format_sumstats #### 
            reformatted <- format_sumstats(
                path = vcf_paths$save_path,
                save_path = save_path,
                log_folder = log_folder,
                force_new = force_new,
                nThread = nThread,
                ...
            )
            end <- Sys.time()
            message(
                "\n", id, " : Done in ",
                round(difftime(end, start, units = "mins"), 3), " minutes."
            )
            return(reformatted)
        }, 
        error = function(e){message(e);e}
        ) 
        return(out)
    }, mc.cores = nThread_acrossIDs) %>% `names<-`(ids)
    
    end_all <- Sys.time()
    message(
        "\nDone with all processing in ",
        round(difftime(end_all, start_all, units = "mins"), 2), " minutes."
    )
    return(ouputs)
}