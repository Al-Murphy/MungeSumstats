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
                           ref_genome,
                           sep = "\t",
                           write_vcf = FALSE,
                           save_format = NULL,
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
    if(isTRUE(tabix_index) |
       isTRUE(write_vcf)) {
      sumstats_dt <- sort_coords(sumstats_dt = sumstats_dt)
    }
    #### Select write format ####
    if (isTRUE(write_vcf)) { 
        tmp_save_path <- gsub("\\.bgz|\\.gz","",save_path)
        #convert to IEU OpenGWAS VCF format (column naming and RSID position)
        if(!is.null(save_format) && 
           tolower(save_format)=="opengwas"){
          #first check genome build - all of openGWAS is GRCh37 currently so 
          #warn user if their data isn't
          gen_build_err <- paste0("Your sumstats has been built on the ",
                                  ref_genome, " reference genome. However, ",
                                  "OpenGWAS is all\n",
                                  "currently built on GRCh37. Use the ",
                                  "convert_ref_genome parameter to liftover to",
                                  " GRCh37 by rerunning format_sumstats.")
          if(toupper(ref_genome)!="GRCH37")
            warning(gen_build_err)
          #necessary cols (https://github.com/MRCIEU/gwas-vcf-specification):
          #NS:NC:ES:SE:LP:AF:AC
          #p-val needs to be -log10 p
          #SNP -> RSID in INFO col
          if("SNP" %in% colnames(sumstats_dt)){
            setnames(sumstats_dt,"SNP","RSID")
          } else{
            stop("SNP/RSID is required for IEU OpenGWAS format VCFs")
          }
          #remove any extra columns
          opengwas_cols <- c("RSID","CHR","BP","A1","A2","P",
                              "FRQ","BETA","SE","N","N_CAS")
          if(any(!(colnames(sumstats_dt) %in% opengwas_cols)))
            sumstats_dt[,colnames(sumstats_dt)[!(colnames(sumstats_dt) %in% 
                                                    opengwas_cols)]:=NULL]
          #### Convert to VRanges and save ####
          vr <- to_vranges(sumstats_dt = sumstats_dt) 
          if("P" %in% names(mcols(vr))){
            #P -> LP
            vr$LP <- -log(vr$P,base=10)
            vr$P <- NULL
          }else{
            stop("P-value is required for IEU OpenGWAS format VCFs")
          }
          #FRQ -> AF
          if("FRQ" %in% names(mcols(vr))){
            vr$AF <- vr$FRQ
            vr$FRQ <- NULL
          }
          #BETA -> ES
          if("BETA" %in% names(mcols(vr))){
            vr$ES <- vr$BETA
            vr$BETA <- NULL
          }else{
            stop("BETA (Effect Size) is required for IEU OpenGWAS format VCFs")
          }
          #SE -> SE
          if(!"SE" %in% names(mcols(vr))){
            se_msg <- paste0("Standard Error (of effect size) is required for",
                              " IEU OpenGWAS format VCFs")
            stop(se_msg)
          }
          #N -> NS
          if("N" %in% names(mcols(vr))){
            vr$NS <- vr$N
            vr$N <- NULL
          }
          #N_CAS -> NC
          if("N_CAS" %in% names(mcols(vr))){
            vr$NC <- vr$N_CAS
            vr$N_CAS <- NULL
          }
          #reorder cols and drop any unnecessary
          all_poss_cols <- c("RSID","NS","NC","ES","SE","LP","AF","AC")
          vr <- vr[,all_poss_cols[all_poss_cols %in% names(mcols(vr))]]
          
          messager("Writing in IEU OpenGWAS VCF format ==> ", save_path)
          VariantAnnotation::writeVcf(
            obj = VariantAnnotation::asVCF(vr,info='RSID'),
            ### Must supply filename without compression suffix
            filename = tmp_save_path,
            index = tabix_index
          )
          
        }else{
          #### Convert to VRanges and save ####
          vr <- to_vranges(sumstats_dt = sumstats_dt)
          messager("Writing in VCF format ==> ", save_path)
          VariantAnnotation::writeVcf(
              obj = vr,
              ### Must supply filename without compression suffix
              filename = tmp_save_path,
              index = tabix_index
          )
        }
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
        if (isTRUE(tabix_index)) {
            tmp_save_path <- gsub(".bgz|.gz", "", save_path)
            messager("Writing in tabular format ==>", tmp_save_path)
            if (tmp_save_path != save_path) {
              messager("Writing uncomressed instead of gzipped to enable index")
            }
        } else {
            tmp_save_path <- save_path
            messager("Writing in tabular format ==>", save_path)

        }
        #### Write to disk #### 
        if (sep == "\t") {
          quotation <- FALSE
          narep <- "NA"
        } else {
          quotation <- "auto"
          narep <- ""
        }

        data.table::fwrite(
            x = sumstats_dt, 
            ### Must be uncompressed so #### 
            file = tmp_save_path,
            sep = sep,
            nThread = nThread,
            na = narep,
            quote = quotation
        )

            save_path <- index_tabular(path = tmp_save_path,
                                       chrom_col = "CHR", 
                                       start_col = "BP", 
                                       verbose = TRUE)
        } 
    }
    if(isTRUE(return_path)) return(save_path)
}