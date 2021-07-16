#' Ensure A1 & A2 are correctly named, if GWAS constructed as 
#' Risk/Nonrisk alleles this will need to be converted to Reference/Alternative. 
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS.
#' @param path Filepath for the summary statistics file to be formatted.
#' @param ref_genome name of the reference genome used for
#'  the GWAS ("GRCh37" or "GRCh38").
#' @param rsids \code{data.table} of snpsById, filtered to SNPs of 
#' interest if loaded already. Or else NULL.
#' @param allele_flip_check Binary Should the allele columns be chacked against 
#' reference genome to infer if flipping is necessary. Default is TRUE.
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
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table set
#' @importFrom data.table setorder
#' @importFrom data.table copy
#' 
#' @examples
#' path <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
#' sumstats_dt <- MungeSumstats::read_sumstats(path = path)
#' sumstats_return <- check_allele_flip(sumstats_dt = sumstats_dt, 
#'                                      path=path, 
#'                                      ref_genome="GRCh37",
#'                                      rsids=NULL,
#'                                      allele_flip_check=TRUE,
#'                                      standardise_headers=TRUE) 
check_allele_flip <-  function(sumstats_dt, 
                               path, 
                               ref_genome, 
                               rsids, 
                               allele_flip_check,
                               keepEA=FALSE,
                               standardise_headers=FALSE){ 
  # GenomicSEM' allele flipping strategy:
  # https://github.com/GenomicSEM/GenomicSEM/blob/fc8f17a817a8022d6900acf41824d27b3676f9c4/R/munge.R#L151
  
  
  message("Checking for alleles to flip.")
  ## Set variables to be used in inplace data.table functions to NULL 
  ## to avoid confusing BiocCheck.
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = A1 = A2 = eff_i =
    i.A1 = i.A2 = ss_A1 = ss_A2 = NULL;
  
  if(standardise_headers){
    sumstats_dt <- standardise_sumstats_column_headers_crossplatform(sumstats_dt = sumstats_dt)[["sumstats_dt"]]
  }
  # If SNP present but no A1/A2 then need to find them
  col_headers <- names(sumstats_dt)
  if(sum(c("A1","A2") %in% col_headers)==2 && allele_flip_check){
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(data.table::copy(sumstats_dt$SNP), ref_genome, NULL)
    }
    #ensure rsids is up-to-date with filtered sumstats_dt
    rsids <- rsids[unique(sumstats_dt$SNP),,nomatch=NULL]
    data.table::setkey(rsids,SNP)
    #Check if reference genome ref allele matches A1 or A2 better in the data
    #If it matches A2 better, flip and flip effect too!
    # join based on SNP as key
    data.table::setorder(sumstats_dt,SNP)
    data.table::setkey(sumstats_dt,SNP)

    rsids[sumstats_dt,ss_A1:=i.A1]
    rsids[sumstats_dt,ss_A2:=i.A2]
    
    ssA1_match <-sum(rsids$ref_allele==rsids$ss_A1)
    ssA2_match <-sum(rsids$ref_allele==rsids$ss_A2)
    
    rsids[, ss_A1:=NULL]
    rsids[, ss_A2:=NULL]
    data.table::setorder(rsids,SNP)
    data.table::setkey(rsids,SNP)
    
    if(ssA2_match>=ssA1_match){#more matches opp dir so flip
      print_msg <- paste0("There are more matches to the reference genome",
                          " found when the alleles are flipped.\nThese will be",
                          " flipped along with the effect column.")
      message(print_msg)
      data.table::setnames(sumstats_dt,"A1","tmp")
      data.table::setnames(sumstats_dt,"A2","A1")
      data.table::setnames(sumstats_dt,"tmp","A2")
      #flip effect column(s) - BETA, OR, z, log_odds, SIGNED_SUMSTAT
      effect_columns <- c("BETA","OR","Z","LOG_ODDS","SIGNED_SUMSTAT")
      effect_columns <- effect_columns[effect_columns %in% names(sumstats_dt)]
      for(eff_i in effect_columns){#set updates quicker for DT
        sumstats_dt[,(eff_i):=get(eff_i)*-1]
      }  
      return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
    }  
  }
  return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
}
