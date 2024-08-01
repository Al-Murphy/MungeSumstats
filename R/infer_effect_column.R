#' Infer if effect relates to a1 or A2 if ambiguously named
#'
#' Three checks are made to infer which allele the effect/frequency information
#' relates to if they are ambiguous (named A1 and A2 or equivalent):
#' 1. Check if ambiguous naming conventions are used (i.e. allele 1 and 2 or 
#' equivalent). If not exit, otherwise continue to next checks. This can be 
#' checked by using the mapping file and splitting A1/A2 mappings by those that 
#' contain 1 or 2 (ambiguous) or doesn't contain 1 or 2 e.g. effect, 
#' tested (unambiguous so fine for MSS to handle as is).
#' 2. Look for effect column/frequency column where the A1/A2 explicitly 
#' mentioned, if found then  we know the direction and should update A1/A2 
#' naming so A2 is the effect column. We can look for such columns by getting 
#' every combination of A1/A2 naming and effect/frq naming.
#' 3. If not found in 2, a final check should be against the reference genome, 
#' whichever of A1 and A2 has more of a match with the reference genome should 
#' be taken as **not** the effect allele. There is an assumption in this but is 
#' still better than guessing the ambiguous allele naming.
#'
#' @inheritParams format_sumstats 
#' @inheritParams compute_nsize
#' @inheritParams get_genome_build
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object
#' @export
#' @importFrom data.table setnames as.data.table := setkey rbindlist data.table
#' @examples
#' sumstats <- MungeSumstats::formatted_example()
#' #for speed, don't run on_ref_genome part of check (on_ref_genome = FALSE)
#' sumstats_dt2<-infer_effect_column(sumstats_dt=sumstats,on_ref_genome = FALSE)
infer_effect_column <-
  function(sumstats_dt,
           dbSNP=155,
           sampled_snps = 10000,
           mapping_file = sumstatsColHeaders,
           nThread = nThread,
           ref_genome = NULL,
           on_ref_genome = TRUE,
           infer_eff_direction = TRUE,
           return_list=TRUE) {
    if(isTRUE(infer_eff_direction)){
      message("Infer Effect Column")
      message("First line of summary statistics file: ")
      msg <- paste0(names(sumstats_dt), split = "\t")
      message(msg)
      #### first make all column headers upper case ####
      column_headers <- names(sumstats_dt) 
      # load synonym mapping - internal data no loading
      # Identify allele mappings which are ambiguous and problematic
      # vs those that are interpretable
      colnames(mapping_file) <- toupper(colnames(mapping_file))
      allele_mapping <- mapping_file[mapping_file$CORRECTED %in% c('A1','A2'),]
      ambig_allele_map <- 
        allele_mapping[grepl('1',allele_mapping$UNCORRECTED)|
                         grepl('2',allele_mapping$UNCORRECTED),]
      unambig_allele_map <- 
        allele_mapping[!(grepl('1',allele_mapping$UNCORRECTED)|
                           grepl('2',allele_mapping$UNCORRECTED)),]
      #as long as the sumstats contains 1 unambiguous allele column MSS will
      #work as expected
      unambig_cols <- intersect(unambig_allele_map$UNCORRECTED,
                                toupper(column_headers))
      ambig_cols <- intersect(ambig_allele_map$UNCORRECTED,
                              toupper(column_headers))
      #if both ambiguous and unambiguous columns found, rename ambiguous ones so
      #they aren't used later by MSS
      #example: 'A1','A2','EFFECT_ALLELE' all present
      if (length(unambig_cols)>0 && length(ambig_cols)>0){
        #find if unambig and ambig relate to the same allele
        #get corrected name for unambig 
        unambig_corrcted <- 
          unique(allele_mapping[allele_mapping$UNCORRECTED %in% unambig_cols,
                                ]$CORRECTED)
        #check if any ambig are to the same allele
        ambig_corrcted <- 
          unique(allele_mapping[allele_mapping$UNCORRECTED %in% ambig_cols,
                                ]$CORRECTED)
        #overlap?
        ambig_corrcted_rnme <- 
          ambig_corrcted[ambig_corrcted %in% unambig_corrcted]
        if (length(ambig_corrcted_rnme)>0){
          message("There are multiple columns relating to the same allele info")
          message("Renaming ambiguous allele columns so they won't be used")
          #get the related ambig naming and change there name so won't be used
          ambig_uncorrcted_rnme <- 
            ambig_allele_map[ambig_allele_map$CORRECTED %in% 
                               ambig_corrcted_rnme,]$UNCORRECTED
          #now rename any matches in sumstats
          chng_nmes <- column_headers[toupper(column_headers) %in% 
                                       ambig_uncorrcted_rnme]
          for(chng_i in chng_nmes){
            data.table::setnames(sumstats_dt, chng_i, 
                                 paste0(chng_i,"_INPUTTED"))
          }
        }
      } else if (length(unambig_cols)==0 && length(ambig_cols)>=2){
        #only continue if no unambiguous columns found but 2 ambig ones are found-
        #less than 2 in total means allele info is missing which MSS can try fill
        #in later
        message("Allele columns are ambiguous, attempting to infer direction")
        #get names for allele marked eff/frq columns
        eff_frq_allele_matches <- get_eff_frq_allele_combns()
        #now look for matches in sumstats
        fnd_allele_indicator <- 
          column_headers[toupper(column_headers) %in% 
                           eff_frq_allele_matches$UNCORRECTED]
        if(length(fnd_allele_indicator)>0){
          message("Found direction from effect/frq column naming")
          #fnd_allele_indicator could be >1 so majority vote
          a1_mtch <- sum(grepl("A1",fnd_allele_indicator))
          a2_mtch <- sum(grepl("A2",fnd_allele_indicator))
          if(a2_mtch>=a1_mtch){
            message("Effect/frq column(s) relate to A2 in the sumstats")
            #this is what MSS expects so no action required
          }else{#a2_mtch<a1_mtch
            message("Effect/frq column(s) relate to A1 in the sumstats")
            #this is the opposite to what MSS expects so switch A1/A2 naming
            #first get corrected names for allele columns then switch
            for (headerI in seq_len(nrow(mapping_file))) {
              un <- mapping_file[headerI, "UNCORRECTED"]
              cr <- mapping_file[headerI, "CORRECTED"]
              #note ambig_cols here not all col headers
              if (un %in% ambig_cols & (!cr %in% column_headers)) {
                data.table::setnames(sumstats_dt, un, cr)
              }
            }
            #now switch
            data.table::setnames(sumstats_dt,"A2","A2_INPUTTED_OLD_")
            data.table::setnames(sumstats_dt,"A1","A2")
            data.table::setnames(sumstats_dt,"A2_INPUTTED_OLD_","A1")
          }
          
        }
        else{
          #didn't find any allele specific
          #last option is to check the reference genome and take the allele with
          #more matches to it as the non-effect allele
          #need to standardise the sumstats for this
          #use get_genome_build function for this
          if(isTRUE(on_ref_genome)){
            switch_req <- get_genome_build(
              sumstats = sumstats_dt,
              standardise_headers = TRUE,
              sampled_snps = sampled_snps,
              mapping_file = mapping_file,
              dbSNP=dbSNP,
              nThread = nThread,
              allele_match_ref = TRUE,
              ref_genome = ref_genome
            )
            if(is.logical(switch_req)){
              message(paste0("Found direction from matching reference genome -",
                             " NOTE this assumes non-effect allele will match ",
                             "the reference genome"))
              if(isTRUE(switch_req)){
                #swap A1 and A2
                #this is the opposite to what MSS expects so switch A1/A2 naming
                #first get corrected names for allele columns then switch
                for (headerI in seq_len(nrow(mapping_file))) {
                  un <- mapping_file[headerI, "UNCORRECTED"]
                  cr <- mapping_file[headerI, "CORRECTED"]
                  #note ambig_cols here not all col headers
                  if (un %in% ambig_cols & (!cr %in% column_headers)) {
                    data.table::setnames(sumstats_dt, un, cr)
                  }
                }
                #now switch
                data.table::setnames(sumstats_dt,"A2","A2_INPUTTED_OLD_")
                data.table::setnames(sumstats_dt,"A1","A2")
                data.table::setnames(sumstats_dt,"A2_INPUTTED_OLD_","A1")
              }
            }
            else{ #couldn't infer
              message(switch_req)
              message("Can't infer allele columns from sumstats")
            }
          }
          else{
            message("Can't infer allele columns from sumstats")
          }
        }
      }
    }  
    #### Return format ####
    if(return_list){
      return(list("sumstats_dt" = sumstats_dt))
    }else {
      return(sumstats_dt)
    }
  }    
