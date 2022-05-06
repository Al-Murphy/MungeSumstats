#' Write sum stats file to disk
#'
#' @param sumstats_dt data table obj of the summary statistics
#' file for the GWAS.
#' @param return_path Return \code{save_path}.
#'  This will have been modified in some cases
#'   (e.g. after compressing and tabix-indexing a
#'    previously un-compressed file).
#' @param save_path_check Ensure path name is valid (given the other arguments) 
#' before writing (default: FALSE). 
#' @inheritParams data.table::fread
#' @inheritParams format_sumstats
#' 
#' @source \href{https://github.com/Bioconductor/VariantAnnotation/issues/35}{
#' VariantAnnotation::writeVcf has some unexpected/silent 
#' file renaming behavior}
#'
#' @returns If \code{return_path=TRUE}, returns \code{save_path}.
#'  Else returns \code{NULL}.
#'
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom VariantAnnotation makeVRangesFromGRanges
#' @examples
#' path <- system.file("extdata", "eduAttainOkbay.txt",
#'     package = "MungeSumstats"
#' )
#' eduAttainOkbay <- read_sumstats(path = path)
#' write_sumstats(
#'     sumstats_dt = eduAttainOkbay,
#'     save_path = tempfile(fileext = ".tsv.gz")
#' )
write_sumstats <- function(sumstats_dt,
                           save_path,
                           sep = "\t",
                           write_vcf = FALSE,
                           tabix_index = FALSE,
                           nThread = 1,
                           return_path = FALSE,
                           save_path_check = FALSE) {
    #### Check save_path ####
    if(isTRUE(save_path_check)){
        check <- check_save_path(save_path = save_path, 
                                 log_folder = tempdir(),
                                 log_folder_ind = FALSE,
                                 tabix_index = tabix_index, 
                                 write_vcf = write_vcf,
                                 verbose = TRUE)
        save_path <- check$save_path
    } 
    #### Sort again just to be sure when tabix-indexing ####
    if(isTRUE(tabix_index) | isTRUE(write_vcf)) {
        sumstats_dt <- sort_coords(sumstats_dt = sumstats_dt)
    }
    #### Select write format ####
    if (isTRUE(write_vcf)) { 
        tmp_save_path <- gsub("\\.bgz|\\.gz","",save_path)
        #### Convert to VRanges and save ####
        vr <- to_vranges(sumstats_dt = sumstats_dt) 
        messager("Writing in VCF format ==> ", save_path)
        VariantAnnotation::writeVcf(
            obj = vr,
            ### Must supply filename without compression suffix
            filename = tmp_save_path,
            index = tabix_index
        )
        #### Compress ####
        ## only compress if this was not already handled by writeVcf(index=T)
        if(isFALSE(tabix_index)){
            if(endsWith(save_path,".bgz")){
                save_path <- Rsamtools::bgzip(tmp_save_path, overwrite=TRUE)
            } else if(endsWith(save_path,".gz")){
                save_path <- R.utils::gzip(tmp_save_path, overwrite=TRUE)
            } else {
                save_path <- tmp_save_path 
            }
        }
    } else { 
        messager("Writing in tabular format ==> ", save_path)
        if(isTRUE(tabix_index)){
            tmp_save_path <- gsub(".bgz|.gz","",save_path)
        } else {
            tmp_save_path <- save_path
        }
        #### Write to disk #### 
        data.table::fwrite(
            x = sumstats_dt, 
            ### Must be uncompressed so #### 
            file = tmp_save_path,
            sep = sep,
            nThread = nThread
        )
        if(isTRUE(tabix_index)){
            save_path <- index_tabular(path = tmp_save_path,
                                       chrom_col = "CHR", 
                                       start_col = "BP", 
                                       verbose = TRUE)
        } 
    }
    if(return_path) return(save_path)
}