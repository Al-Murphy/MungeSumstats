#' Ensure all SNPs are on the reference genome
#'
#' @inheritParams format_sumstats  
#' @param log_files list of log file locations
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics data table object
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#'   \item \code{log_files}: log file list
#' }
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table copy
#' @importFrom data.table rbindlist
#' @importFrom data.table setcolorder
check_on_ref_genome <- 
  function(sumstats_dt, path, ref_genome, on_ref_genome, rsids,imputation_ind,
           log_folder_ind,check_save_out,tabix_index, nThread,log_files){
  CHR = SNP = IMPUTATION_SNP = miss_rs_chr_bp = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_dt)
  if("SNP" %in% col_headers && !isFALSE(on_ref_genome)){
    message("Ensuring all SNPs are on the reference genome.")
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(snps = data.table::copy(sumstats_dt$SNP), 
                                    ref_genome = ref_genome)
    }
    #ensure rsids is up-to-date with filtered sumstats_dt
    rsids <- rsids[unique(sumstats_dt$SNP),,nomatch=NULL]
    data.table::setkey(rsids,SNP)
    num_bad_ids <- length(sumstats_dt$SNP) - length(rsids$SNP)
    #check for SNPs not on ref genome
    if(num_bad_ids>0){
      msg <- paste0(formatC(num_bad_ids,big.mark = ","), 
                    " SNPs are not on the reference genome.",
                    " These will be corrected from the reference genome.")
      message(msg)
      # join using SNP
      data.table::setkey(sumstats_dt,SNP)
      #if the dataset has CHR & BP, try impute the correct ones
      if(sum(c("CHR","BP") %in% col_headers)==2){
        bad_snp <- sumstats_dt[!rsids$SNP,]
        #remove snp column and pass to function to impute snp
        bad_snp <- bad_snp[,SNP:=NULL]
        corrected_snp <- 
          check_no_snp(sumstats_dt=bad_snp, path=tempfile(), 
                       ref_genome=ref_genome, imputation_ind=imputation_ind,
                       log_folder_ind=log_folder_ind,
                       check_save_out=check_save_out,tabix_index = tabix_index,
                       nThread = nThread,log_files = log_files,verbose=FALSE)
        log_files <- corrected_snp$log_files
        corrected_snp <- corrected_snp$sumstats_dt 
        #make sure columns in correct order
        data.table::setcolorder(corrected_snp,names(sumstats_dt))
        #remove rows missing from the reference genome and combine
        #If IMPUTATION column added add it to other DT
        if(imputation_ind && !"IMPUTATION_SNP" %in% names(sumstats_dt))
          sumstats_dt[,IMPUTATION_SNP:=NA]
        sumstats_dt <- 
          data.table::rbindlist(list(sumstats_dt[rsids$SNP,],corrected_snp))
      }
      else{
        msg <- paste0(formatC(num_bad_ids,big.mark = ","), 
                      " SNPs are not on the reference genome. ",
                      "These will be removed")
        message(msg)
        
        #remove rows missing from the reference genome
        #If user wants log, save it to there
        if(log_folder_ind){
          name <- "snp_not_on_ref_gen"
          name <- get_unique_name_log_file(name=name,log_files=log_files)
          write_sumstats(sumstats_dt = sumstats_dt[!rsids$SNP,],
                         save_path=
                           paste0(check_save_out$log_folder,
                                  "/",name,
                                  check_save_out$extension),
                         sep=check_save_out$sep,
                         tabix_index = tabix_index,
                         nThread = nThread)
          log_files[[name]] <- 
            paste0(check_save_out$log_folder,"/",name,check_save_out$extension)
        } 
        sumstats_dt <- sumstats_dt[rsids$SNP,]
      }
    }
  }
  return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids,"log_files"=log_files))
}
