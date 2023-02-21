#' Ensure that the BAse-apir column values are all within the rnage for the 
#' chromosome
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
#' @importFrom GenomeInfoDb seqinfo
check_bp_range <- function(sumstats_dt, path, ref_genome,log_folder_ind, 
                         imputation_ind,check_save_out, tabix_index, nThread, 
                         log_files) {
    chr_max<- CHR <- seqlengths <- .SD <- NULL
    col_headers <- names(sumstats_dt)
    if ("BP" %in% col_headers) {
        message("Checking for incorrect base-pair positions")
        #load chromosome lengths
        if (toupper(ref_genome) == "GRCH37") {
          genome_info <- 
            seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)
        } else { # =="GRCH38"
          genome_info <- 
            seqinfo(BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38)
        }
        #mke data.table
        genome_info <-
          cbind(data.table::data.table("CHR"=seqnames(genome_info)),
                data.table::as.data.table(genome_info))
        #convwert CHR to char in sumstats (if it is an int)
        sumstats_dt[,CHR:=as.character(CHR)]
        #add chr len to sumstats
        sumstats_dt[genome_info,chr_max:=seqlengths,on="CHR"]
        #now check if BP value is too large/less than 0
        #s.na(chr_max) in case chr formatting not found
        bad_bp <- sumstats_dt[!is.na(chr_max) & (is.na(BP)|BP<=0 |BP>chr_max),]
        if (nrow(bad_bp)>0){
          message(nrow(bad_bp)," SNPs have been removed as their BP column is ",
                  "not in the range of 1 to the length of the chromosome")
          if (log_folder_ind){
            name <- "bad_bp"
            name <- get_unique_name_log_file(name = name,
                                             log_files = log_files)
            write_sumstats(
              sumstats_dt =
                bad_bp,
              save_path =
                paste0(
                  check_save_out$log_folder,
                  "/", name,
                  check_save_out$extension
                ),
              sep = check_save_out$sep,
              tabix_index = tabix_index,
              nThread = nThread
            )
            log_files[[name]] <-
              paste0(check_save_out$log_folder, "/",
                     name, check_save_out$extension)
          }
          
          # update actual
          sumstats_dt <- sumstats_dt[!is.na(BP) & BP>0 & 
                                       (BP<=chr_max | is.na(chr_max)),]
        }
        #check is doen remove column
        sumstats_dt[,chr_max:=NULL]
    }
    
    return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
}