#' Ensure A1 & A2 are correctly named, if GWAS SNP constructed as 
#' Alternative/Reference or Risk/Nonrisk alleles these SNPs will need to be 
#' converted to Reference/Alternative or Nonrisk/Risk. Here nonrisk is defined
#' as what's on the reference genome (this may not always be the case). 
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @param path Filepath for the summary statistics file to be formatted.
#' @param ref_genome name of the reference genome used for
#'  the GWAS ("GRCh37" or "GRCh38").
#' @param rsids \code{data.table} of snpsById, filtered to SNPs of 
#' interest if loaded already. Or else NULL.
#' @param allele_flip_check Binary Should the allele columns be checked against 
#' reference genome to infer if flipping is necessary. Default is TRUE.
#' @param allele_flip_drop Binary Should the SNPs for which neither their A1 or 
#' A2 base pair values match a reference genome be dropped. Default is TRUE.
#' @param keepEA Add a column "EA" to record which allele was the effect allele 
#' in the original summary statistics.
#' @param standardise_headers Run \code{standardise_sumstats_column_headers_crossplatform} first.  
#' 
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics \code{data.table} object.
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest if loaded already. Or else NULL.
#' }
#' @keywords internal
#' @import R.utils
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table set
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_allele_flip <-  function(sumstats_dt, 
                               path, 
                               ref_genome, 
                               rsids, 
                               allele_flip_check,
                               allele_flip_drop,
                               allele_flip_z,
                               keepEA=FALSE,
                               standardise_headers=FALSE){ 
  # GenomicSEM' allele flipping strategy:
  # https://github.com/GenomicSEM/GenomicSEM/blob/fc8f17a817a8022d6900acf41824d27b3676f9c4/R/munge.R#L151
  
  # #example
  # path <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
  # sumstats_dt <- MungeSumstats::read_sumstats(path = path)
  # sumstats_return <- check_allele_flip(sumstats_dt = sumstats_dt,
  #                                      path=path,
  #                                      ref_genome="GRCh37",
  #                                      rsids=NULL,
  #                                      allele_flip_check=TRUE,
  #                                      standardise_headers=TRUE)
  ## Set variables to be used in in place data.table functions to NULL 
  ## to avoid confusing BiocCheck.
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = A1 = A2 = eff_i =
    i.A1 = i.A2 = ss_A1 = ss_A2 = i.ref_allele = ref_gen_allele = match_type = 
    tmp = NULL;
  if(standardise_headers){
    sumstats_dt <- standardise_sumstats_column_headers_crossplatform(sumstats_dt = sumstats_dt)[["sumstats_dt"]]
  }
  # If SNP present but no A1/A2 then need to find them
  col_headers <- names(sumstats_dt)
  if(sum(c("A1","A2") %in% col_headers)==2 && allele_flip_check){
    message("Checking for correct direction of A1 (reference) and A2 (alternative allele).")
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- 
        load_ref_genome_data(data.table::copy(sumstats_dt$SNP),ref_genome, NULL)
    }
    #ensure rsids is up-to-date with filtered sumstats_dt
    rsids <- rsids[unique(sumstats_dt$SNP),,nomatch=NULL]
    data.table::setkey(rsids,SNP)
    data.table::setkey(sumstats_dt,SNP)
    #Append reference genome to data and check matches to both A1 or A2 
    # For each SNP, if ref genome matches A1, leave it (TRUE)
    # if ref genome matches A2, flip it (FALSE)
    # if ref genome doesn't match A1 or A2, leave it (TRUE)
    sumstats_dt[rsids,ref_gen_allele:=i.ref_allele]
    sumstats_dt[is.na(ref_gen_allele),match_type:=TRUE]
    sumstats_dt[A1==ref_gen_allele,match_type:=TRUE]
    sumstats_dt[A2==ref_gen_allele,match_type:=FALSE]
    #drop cases that don't match either
    if(allele_flip_drop && 
        nrow(sumstats_dt[A1!=ref_gen_allele & A2!=ref_gen_allele,])>0){
      print_msg0 <- 
        paste0("There are ",
               nrow(sumstats_dt[A1!=ref_gen_allele & A2!=ref_gen_allele,]),
               " SNPs where neither A1 nor A2 match the reference genome.",
               "\nThese will be removed.")
      message(print_msg0)
      sumstats_dt <- sumstats_dt[!(A1!=ref_gen_allele & A2!=ref_gen_allele),]
    }
    else{
      sumstats_dt[A1!=ref_gen_allele & A2!=ref_gen_allele,match_type:=TRUE]  
    }
    #continue if flipping necessary
    if(nrow(sumstats_dt[match_type==FALSE,])>0){
      print_msg <- paste0("There are ",nrow(sumstats_dt[match_type==FALSE,]),
                          " SNPs where A1 doesn't match the reference genome.",
                          "\nThese will be flipped with their effect columns.")
      message(print_msg)
      #swap A1 and A2 for those SNPs needing to be flipped
      sumstats_dt[match_type==FALSE,tmp:=A2]
      sumstats_dt[match_type==FALSE,A2:=A1]
      sumstats_dt[match_type==FALSE,A1:=tmp]
      sumstats_dt[,tmp:=NULL]
      
      #flip effect column(s) - BETA, OR, z, log_odds, SIGNED_SUMSTAT
      effect_columns <- c("BETA","OR","LOG_ODDS","SIGNED_SUMSTAT")
      if(allele_flip_z)
        effect_columns <- c("BETA","OR","Z","LOG_ODDS","SIGNED_SUMSTAT")
      effect_columns <- effect_columns[effect_columns %in% names(sumstats_dt)]
      for(eff_i in effect_columns){#set updates quicker for DT
        #conversion done in case, VCF beta column may not be numeric
        sumstats_dt[match_type==FALSE,(eff_i):=as.numeric(get(eff_i))*-1]
      }  
    }
    #remove extra created columns and return
    sumstats_dt[,ref_gen_allele:=NULL]
    sumstats_dt[,match_type:=NULL]

    return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
  }  
  return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
}
