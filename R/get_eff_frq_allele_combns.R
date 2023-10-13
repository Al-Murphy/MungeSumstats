#' Get combinations of uncorrected allele and effect (and frq) columns
#'
#' @inheritParams format_sumstats 
#' @inheritParams compute_nsize
#' @param eff_frq_cols Corrected effect or frequency column names found in a 
#' sumstats. Default of BETA, OR, LOG_ODDS, SIGNED_SUMSTAT, Z and FRQ.
#' @return datatable containing uncorrected and corrected combinations
#' @importFrom data.table setnames as.data.table := setkey rbindlist data.table
get_eff_frq_allele_combns <-
  function(mapping_file = sumstatsColHeaders, 
           eff_frq_cols = c("BETA", "OR", "LOG_ODDS", "SIGNED_SUMSTAT","Z",
                            "FRQ")) {
    ### Add this to avoid confusing BiocCheck
    CORRECTED <- UNCORRECTED <- Var1 <- Var2 <- NULL 
    colnames(mapping_file) <- toupper(colnames(mapping_file))
    #get allele associated effect/FRQ columns
    #get all combinations with allele columns
    eff_frq_cols_uncorrc <- 
      mapping_file[mapping_file$CORRECTED %in% eff_frq_cols,]$UNCORRECTED
    #join with all allele cols
    allele_uncorrc <- 
      mapping_file[mapping_file$CORRECTED %in% c('A1','A2'),]$UNCORRECTED
    #get combinations
    eff_frq_allele_dt <- 
      data.table::as.data.table(expand.grid(eff_frq_cols_uncorrc, 
                                            allele_uncorrc))
    mapping_file_dt <- data.table::as.data.table(mapping_file)
    data.table::setkey(mapping_file_dt,"UNCORRECTED")
    data.table::setkey(eff_frq_allele_dt,"Var1")
    #add corrected
    eff_frq_allele_dt[mapping_file,CORRECTED:=CORRECTED,]
    #now loop through every joining character and join with eff both before
    #and after
    joining_char <- c("","_",".","-"," ")
    all_combns <- vector(mode="list",length = length(joining_char)*2)
    counter <- 1
    for(join_i in joining_char){
      eff_frq_allele_dt_i <- copy(eff_frq_allele_dt)  
      eff_frq_allele_dt_i[,UNCORRECTED:=paste0(Var1,join_i,Var2)]
      all_combns[[counter]] <- 
        eff_frq_allele_dt_i[,c("UNCORRECTED","CORRECTED")]
      counter <- counter+1
      #same for Var 2 in front
      eff_frq_allele_dt_i <- copy(eff_frq_allele_dt)  
      eff_frq_allele_dt_i[,UNCORRECTED:=paste0(Var2,join_i,Var1)]
      all_combns[[counter]] <- 
        eff_frq_allele_dt_i[,c("UNCORRECTED","CORRECTED")]
      counter <- counter+1
    }
    #join all together
    eff_frq_allele_matches <- data.table::rbindlist(all_combns)
    #finally add some custom ones 
    custom_adds <- data.table::data.table("UNCORRECTED" = 
                                            c("BETA1", "BETA2","AF1","AF2",
                                              "FREQ.A1.1000G.EUR",
                                              "FREQ.A2.1000G.EUR",
                                              "FREQ.A1.ESP.EUR",
                                              "FREQ.A2.ESP.EUR",
                                              "FREQ.ALLELE1.HAPMAPCEU",
                                              "FREQ.ALLELE2.HAPMAPCEU",
                                              "FREQ1","FREQ2", 
                                              "FREQ1.HAPMAP","FREQ2.HAPMAP"),
                                          "CORRECTED" = 
                                            c("BETA", "BETA","FRQ","FRQ",
                                              "FRQ",
                                              "FRQ",
                                              "FRQ",
                                              "FRQ",
                                              "FRQ",
                                              "FRQ",
                                              "FRQ","FRQ",
                                              "FRQ","FRQ"))
    eff_frq_allele_matches <- data.table::rbindlist(list(
      eff_frq_allele_matches,custom_adds))
    
    return(eff_frq_allele_matches)
  }
