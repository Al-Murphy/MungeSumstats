#' Ensure A1 & A2 are correctly named, if GWAS constructed as Risk/nonrisk alleles this will need to be converted
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param rsids datatable of snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' @param allele_flip_check Binary Should the allele columns be chacked against reference genome to infer if flipping is necessary. Default is TRUE
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics data table object
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' }
#' @keywords internal
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table set
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_allele_flip <- 
  function(sumstats_dt, path, ref_genome, rsids, allele_flip_check){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = A1 = A2 = eff_i =
    i.A1 = i.A2 = ss_A1 = ss_A2 = NULL
  # If SNP present but no A1/A2 then need to find them
  col_headers <- names(sumstats_dt)
  if(sum(c("A1","A2") %in% col_headers)==2 & allele_flip_check){
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_dt$SNP), ref_genome, NULL)
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
