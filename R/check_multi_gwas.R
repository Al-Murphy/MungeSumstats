#' Ensure that only one model in GWAS sumstats or only one trait tested
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param analysis_trait If multiple traits were studied, name of the trait for analysis from the GWAS. Default is NULL
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
check_multi_gwas <- function(sumstats_dt, path, analysis_trait){
  sumstatsColHeaders <- MungeSumstats:::sumstatsColHeaders
  col_headers <- names(sumstats_dt)
  #load synonym mapping - internal data no loading
  #Check for multiple p values and other columns, if found choose 1 to continue
  #Use P value to find - first load all p-value synonyms
  p_syn <- sumstatsColHeaders[sumstatsColHeaders$Corrected=="P", 1]
  # So P_value columns found have already been corrected to "P" so if you find
  # any in p_syn in colheaders then we found multiple
  # NOTE THIS ASSUMES AN '_' or '.' SEPARATES P VALUE SYN FROM TRAIT/MODEL NAME
  # AND COLUMN NAME STARTS WITH P VAULE SYM
  # Sort p_syn so longest ones checked first, once 1 found stop
  # so don't get wrong subset with shorter ones
  p_syn <- p_syn[order(nchar(p_syn), p_syn, decreasing = TRUE)]
  for(pattern_i in p_syn){
    p_val_multi <- grepl(paste0("^",pattern_i,"_","|",
                                  "^",pattern_i,"[.]"),col_headers)
    if(any(p_val_multi))
      break
  }
  # So note we could have 0 (do nothing), 1 (just update names),
  # 2+ (remove all but one) traits
  if (any(p_val_multi)) {#multiple traits
    #first get trait name by removing specific p_syn
    traits <- gsub(paste0("^",pattern_i,"_","|",
                            "^",pattern_i,"[.]"),"",col_headers[p_val_multi])
    if(length(traits)>1){
      if(!is.null(analysis_trait))
        analysis_trait <- as.character(analysis_trait)#ensure character
      msg <- paste0("WARNING: Multiple traits found in sumstats file only one",
                    " of which can be analysed: \n",
                    paste(traits, collapse=', ' ))
      message(msg)
      stp_msg <- paste0("Inputted analysis_trait not one of these trait. Pleas",
                        "e input one of the traits above.")
      if(is.null(analysis_trait) ||
         !toupper(analysis_trait) %in% toupper(traits))
        stop(stp_msg)
      #get right case for choice
      chosen_trait <- traits[toupper(traits) %in% toupper(analysis_trait)]
    }
    else{#just one trait
      chosen_trait <- traits
    }
    #just 1 trait from start or 1 chosen by user, remove the trait from col name
    #then get the correct synonyms
    #NOTE ASSUMES TRAIT NAME STARTS WITH '_' or '.' & IS AT THE END OF COLUMN
    chnge_header_names <- col_headers[grepl(paste0("_",chosen_trait,"$","|",
                                                    "[.]",chosen_trait,"$"),
                                              col_headers)]
    
    new_names <- gsub(paste0("_",chosen_trait,"$","|",
                                ".",chosen_trait,"$"),"",chnge_header_names)
    data.table::setnames(sumstats_dt,chnge_header_names,new_names)
    #RE-standardise headers for all OS
    sumstats_return <-
      standardise_sumstats_column_headers_crossplatform(sumstats_dt, path)

    return(sumstats_return)
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
