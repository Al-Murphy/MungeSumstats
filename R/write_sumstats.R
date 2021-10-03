#' Write sum stats file to disk
#'
#' @param sumstats_dt data table obj of the summary statistics
#' file for the GWAS.
#' @param return_path Return \code{save_path}.
#'  This will have been modified in some cases
#'   (e.g. after compressing and tabix-indexing a
#'    previously un-compressed file).
#' @inheritParams data.table::fread
#' @inheritParams format_sumstats
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
                           return_path = FALSE) {
    #### Make sure the directory actually exists
    if (is.character(save_path)) {
        dir.create(dirname(save_path),
            showWarnings = FALSE,
            recursive = TRUE
        )
    }
    if (write_vcf) {
        vr <- to_vranges(sumstats_dt = sumstats_dt)
        if (tabix_index) {
            suffixes <- supported_suffixes(
                tabular = FALSE,
                tabular_compressed = FALSE
            )
            message("Writing in VCF format ==> ", gsub(
                paste(suffixes, collapse = "|"),
                ".vcf.bgz", save_path
            ))
            message("Compressing with bgzip and indexing with tabix.")
        } else {
            message("Writing in VCF format ==> ", save_path)
        }
        VariantAnnotation::writeVcf(
            obj = vr,
            filename = save_path,
            index = tabix_index
        )
    } else {
        message("Writing in tabular format ==> ", save_path)
        data.table::fwrite(
            x = sumstats_dt,
            file = save_path,
            sep = sep,
            nThread = nThread
        )
        if(tabix_index){
            save_path <- index_tabular(path = save_path,
                                       chrom_col = "CHR", 
                                       start_col = "BP", 
                                       verbose = TRUE)
        }
    }
    if(return_path) return(save_path)
}
