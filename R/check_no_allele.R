#' Ensure that A1 & A2 are present, if not can find it with SNP and other allele
#' 
#' More care needs to be taken if one of A1/A2 is present, before imputing the 
#' other allele flipping needs to be checked 
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param rsids datatable of snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics data table object
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' }
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table setcolorder
#' @importFrom data.table setorder
#' @importFrom data.table copy
#' @importFrom Biostrings IUPAC_CODE_MAP 
check_no_allele <- function(sumstats_dt, path, ref_genome, rsids){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = A1 = A2 = 
    i.ref_allele = i.alt_alleles = alt_alleles = ss_A = alleles_as_ambig = NULL
  # If SNP present but no A1/A2 then need to find them
  col_headers <- names(sumstats_dt)
  if(sum(c("A1","A2") %in% col_headers)<=1 & sum("SNP" %in% col_headers)==1){
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_dt$SNP), ref_genome,
                                      "A1 or A2 allele information")
    }
    else{
      print_msg <- paste0("There is no A1 or A2 allele information column",
                          " found within the data. It must be inferred from ",
                          " other column information.")
      message(print_msg)
    }
    
    # join based on SNP as key
    data.table::setkey(sumstats_dt,SNP)
    #if one allele in dataset join other
    if(sum(c("A1","A2") %in% col_headers)==1){ #one allele missing
      # First check if reference genome ref/alt allele matches allele better 
      # If it matches allele not equal to data allele's name better, 
      # flip and flip effect too! e.g. If have A1 but matches A2 in ref flip
      #get col name in data
      ssA <- c("A1","A2")[c("A1","A2") %in% col_headers]
      # join based on SNP as key
      data.table::setorder(sumstats_dt,SNP)
      data.table::setkey(sumstats_dt,SNP)
      
      rsids[sumstats_dt,ss_A:=get(ssA)]
      #exclude bi/tri allelic from check
      #get chars for SNPs not bi/tri allelic or strand ambig from IUPAC_CODE_MAP
      nonambig_IUPAC_CODE_MAP <- 
        names(Biostrings::IUPAC_CODE_MAP[nchar(Biostrings::IUPAC_CODE_MAP)<3])
      A1_match <-sum(rsids[alleles_as_ambig %in% 
                             nonambig_IUPAC_CODE_MAP]$ref_allele==rsids$ss_A)
      #convert to char as currently a list for bi/tri allelic
      A2_match <-sum(as.character(rsids[alleles_as_ambig %in% 
                             nonambig_IUPAC_CODE_MAP]$alt_allele)==rsids$ss_A)
      match_scores <- c("A1"=A1_match,"A2"=A2_match)
      more_matched_allele <- names(which.max(match_scores))
      rsids[, ss_A:=NULL]
      data.table::setorder(rsids,SNP)
      data.table::setkey(rsids,SNP)
      
      if(more_matched_allele!=ssA){
        print_msg <- paste0("There are more matches to the reference genome",
                            " found when the supplied allele is flipped. This ",
                            "will be flipped along with the effect column.")
        message(print_msg)
        #Update name to other allele
        data.table::setnames(sumstats_dt,
                             names(match_scores)[!names(match_scores) %in% 
                                                   more_matched_allele],
                             more_matched_allele)
        #flip effect column(s) - BETA, OR, z, log_odds, SIGNED_SUMSTAT
        effect_columns <- c("BETA","OR","Z","LOG_ODDS","SIGNED_SUMSTAT")
        effect_columns <- 
          effect_columns[effect_columns %in% names(sumstats_dt)]
        for(eff_i in effect_columns){#set updates quicker for DT
          data.table::set(sumstats_dt, i=NULL, j=eff_i, 
                          value=sumstats_dt[[eff_i]]*-1)
        }  
      }
      
      #Now that we have the correct allele name, add other allele from ref gen
      if("A1" %in% names(sumstats_dt)){
        message("Deriving A2 from reference genome")
        msg <- paste0("WARNING: Inferring the alternative allele (A2) from the",
                      " reference genome. In some instances, there are more ",
                      "than one\nalternative allele. Arbitrarily, only the ",
                      "first will be kept. See column `alt_alleles` in your ",
                      "returned sumstats file\nfor all alternative alleles.")
        message(msg)
        sumstats_dt[rsids,alt_alleles:=i.alt_alleles]
        #just take first A2 value arbitrarily
        sumstats_dt[,A2:= as.character(lapply(alt_alleles, function(x) x[1]))]
        #collapse alt_alleles into character type sep by columns
        sumstats_dt[,alt_alleles:= 
                        as.character(lapply(alt_alleles, 
                                            function(x) paste0(x,
                                                               collapse= ",")))]
      }
      else{ #A2 in input
        message("Deriving A1 from reference genome")
        sumstats_dt[rsids,A1:=i.ref_allele]
      }
    }
    else{ #get both A1, A2 from ref genome - choose an A2 value where multiple
      message("Deriving both A1 and A2 from reference genome")
      sumstats_dt[rsids,A1:=i.ref_allele]
      msg <- paste0("WARNING: Inferring the alternative allele (A2) from the ",
                    "reference genome. In some instances, there are more ",
                    "than one\nalternative allele. Arbitrarily, only the ",
                    "first will be kept. See column `alt_alleles` in your ",
                    "returned sumstats file\nfor all alternative alleles.")
      message(msg)
      sumstats_dt[rsids,alt_alleles:=i.alt_alleles]
      #just take first A2 value arbitrarily
      sumstats_dt[,A2:= as.character(lapply(alt_alleles, function(x) x[1]))]
      #collapse alt_alleles into character type sep by columns
      sumstats_dt[,alt_alleles:= 
                      as.character(lapply(alt_alleles, 
                                          function(x) paste0(x,
                                                             collapse = ",")))]
    }
    #remove rows where A1/A2 couldn't be found
    sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt),]
    #move SNP, CHR, BP, A1 and A2 to start
    other_cols <-
      names(sumstats_dt)[!names(sumstats_dt) %in% 
                             c("SNP","CHR","BP", "A1", "A2")]
    data.table::setcolorder(sumstats_dt, 
                              c("SNP","CHR","BP", "A1", "A2", other_cols))
    
    return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
  }
}
