#' Ensure A1 & A2 are correctly named, if GWAS constructed as Risk/nonrisk alleles this will need to be converted
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param rsids datatable of snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' @param allele_flip_check Binary Should the allele columns be chacked against reference genome to infer if flipping is necessary. Default is TRUE
#' @return The modified sumstats_file
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table set
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_allele_flip <- 
  function(sumstats_file, path, ref_genome, rsids, allele_flip_check){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = A1 = A2 = 
    i.A1 = i.A2 = ss_A1 = ss_A2 = NULL
  # If SNP present but no A1/A2 then need to find them
  col_headers <- names(sumstats_file)
  if(sum(c("A1","A2") %in% col_headers)==2 & allele_flip_check){
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_file$SNP), ref_genome, NULL)
      #Save to parent environment so don't have to load again
      assign("rsids", rsids, envir = parent.frame())
    }
    
    #Check if reference genome ref allele matches A1 or A2 better in the data
    #If it matches A2 better, flip and flip effect too!
    # join based on SNP as key
    data.table::setorder(sumstats_file,SNP)
    data.table::setkey(sumstats_file,SNP)
    
    rsids[sumstats_file,ss_A1:=i.A1]
    rsids[sumstats_file,ss_A2:=i.A2]
    
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
      data.table::setnames(sumstats_file,"A1","tmp")
      data.table::setnames(sumstats_file,"A2","A1")
      data.table::setnames(sumstats_file,"tmp","A2")
      #flip effect column(s) - BETA, OR, z, log_odds, SIGNED_SUMSTAT
      effect_columns <- c("BETA","OR","Z","LOG_ODDS","SIGNED_SUMSTAT")
      effect_columns <- effect_columns[effect_columns %in% names(sumstats_file)]
      for(eff_i in effect_columns){#set updates quicker for DT
        data.table::set(sumstats_file, i=NULL, j=eff_i, 
                          value=sumstats_file[[eff_i]]*-1)
      }  
      return(sumstats_file)
    }  
  }
  return(sumstats_file)
}
