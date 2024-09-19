#' Infer if effect relates to a1 or A2 if ambiguously named
#'
#' Three checks are made to infer which allele the effect/frequency information
#' relates to if they are ambiguous (named A0, A1 and A2 or equivalent):
#' 1. Check if ambiguous naming conventions are used (i.e. allele 0, 1 and 2 or 
#' equivalent). If not exit, otherwise continue to next checks. This can be 
#' checked by using the mapping file and splitting A1/A2 mappings by those that 
#' contain 0, 1 or 2 (ambiguous) or doesn't contain 0, 1 or 2 e.g. effect, 
#' tested (unambiguous so fine for MSS to handle as is).
#' 2. Look for effect column/frequency column where the A0/A1/A2 explicitly 
#' mentioned, if found then  we know the direction and should update A0/A1/A2 
#' naming so A2 is the effect column. We can look for such columns by getting 
#' every combination of A0/A1/A2 naming and effect/frq naming.
#' 3. If not found in 2, a final check should be against the reference genome, 
#' whichever of A0, A1 and A2 has more of a match with the reference genome 
#' should be taken as **not** the effect allele. There is an assumption in this 
#' but is still better than guessing the ambiguous allele naming.
#' 
#' Also, if eff_on_minor_alleles=TRUE, check 3 will be used in all cases. 
#' However, This assumes that the effects are majoritively measured on the 
#' minor alleles and should be used with caution as this is an assumption that 
#' won't be appropriate in all cases. However, the benefit is that if we know 
#' the majority of SNPs have their effects based on the minor alleles, we can 
#' catch cases where the allele columns have been mislabelled. IF 
#' eff_on_minor_alleles=TRUE, checks 1 and 2 will be skipped.
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
           eff_on_minor_alleles = FALSE,
           return_list=TRUE) {
    if(isTRUE(infer_eff_direction)||isTRUE(eff_on_minor_alleles)){
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
      #A* is A0, A* used since usually if A0/A1 used, meaning of A1 
      #becomes eff allele and A0 is non-eff so need to flip later
      allele_mapping <- mapping_file[mapping_file$CORRECTED %in% c('A1','A2',
                                                                   'A*'),]
      ambig_allele_map <- 
        allele_mapping[grepl('1',allele_mapping$UNCORRECTED)|
                         grepl('2',allele_mapping$UNCORRECTED)|
                         grepl('0',allele_mapping$UNCORRECTED),]
      unambig_allele_map <- 
        allele_mapping[!(grepl('1',allele_mapping$UNCORRECTED)|
                           grepl('2',allele_mapping$UNCORRECTED)|
                           grepl('0',allele_mapping$UNCORRECTED)),]
      #as long as the sumstats contains 1 unambiguous allele column MSS will
      #work as expected
      unambig_cols <- intersect(unambig_allele_map$UNCORRECTED,
                                toupper(column_headers))
      ambig_cols <- intersect(ambig_allele_map$UNCORRECTED,
                              toupper(column_headers))
      #if both ambiguous and unambiguous columns found, rename ambiguous ones
      #so they aren't used later by MSS
      #example: 'A1','A2','EFFECT_ALLELE' all present
      if (length(unambig_cols)>0 && length(ambig_cols)>0 && 
          isFALSE(eff_on_minor_alleles)){
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
      } else if ((length(unambig_cols)==0 && length(ambig_cols)>=2) || 
                 isTRUE(eff_on_minor_alleles)){
        #first case for ambig allelee where user didn't set eff_on_minor_alleles
        if ((length(unambig_cols)==0 && length(ambig_cols)>=2) && 
            isFALSE(eff_on_minor_alleles)){
          #only continue if no unambiguous columns found but 2 ambig ones are 
          #found- less than 2 in total means allele info is missing which MSS 
          #can try fill in later
          message("Allele columns are ambiguous, attempting to infer direction")
          #get names for allele marked eff/frq columns
          eff_frq_allele_matches <- get_eff_frq_allele_combns()
          #now look for matches in sumstats
          fnd_allele_indicator <- 
            column_headers[toupper(column_headers) %in% 
                             eff_frq_allele_matches$UNCORRECTED]
        } else{
          #for eff_on_minor_alleles = TRUE - 
          #force length(fnd_allele_indicator)>0 to return FALSE
          fnd_allele_indicator<-c()
        }
        if(length(fnd_allele_indicator)>0){
          message("Found direction from effect/frq column naming")
          #fnd_allele_indicator could be >1 so majority vote
          a1_mtch <- sum(grepl("A1",fnd_allele_indicator))
          a2_mtch <- sum(grepl("A2",fnd_allele_indicator))
          a0_mtch <- sum(grepl("A0",fnd_allele_indicator))
          #need to also check if allele 0 & allele 1 found or more normal case
          #of allele 1 & allele 2, as this flips which is interp as the eff 
          #allele by MSS
          samp_dt <- copy(sumstats_dt[1:10])
          samp_dt <- 
            standardise_sumstats_column_headers_crossplatform(samp_dt,
                                                              mapping_file=
                                                                mapping_file,
                                                              convert_A0=FALSE,
                                                              return_list=FALSE)
          formatted_col_headers <- names(samp_dt)
          #check if A0,A1 and A2 present
          a1_found <- "A1" %in% formatted_col_headers
          a2_found <- "A2" %in% formatted_col_headers
          a0_found <- "A*" %in% formatted_col_headers
          #if A2 found at all, it is eff col normally in MSS
          if(a2_mtch>=a1_mtch && a2_found){
            message("Effect/frq column(s) relate to A2 in the sumstat")
            #this is what MSS expects so no action required
          }else if(a2_mtch<a1_mtch && a2_found){
            message("Effect/frq column(s) relate to A1 in the sumstat")
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
          }else if(a1_mtch>=a0_mtch&& a0_found){
            message("Effect/frq column(s) relate to A1 where A0 in the sumstat")
            #this is what MSS expects so no action required
          }else if(a1_mtch<a0_mtch&& a0_found){
            message("Effect/frq column(s) relate to A0 in the sumstat")
            #this is the opposite to what MSS expects so switch A1/A0 naming
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
            data.table::setnames(sumstats_dt,"A*","A0_INPUTTED_OLD_")
            #don't flip A1 as the next step of the mapping is to flip it anyway
            data.table::setnames(sumstats_dt,"A0_INPUTTED_OLD_","A2")
          }else{#unknown
            print("ERROR: Unknown condition, not changing allele mapping.")
          }
          
        }
        else{
          #Either - didn't find any allele specific or set eff_on_minor_alleles
          #last option is to check the reference genome and take the allele with
          #more matches to it as the non-effect allele
          #need to standardise the sumstats for this
          #use get_genome_build function for this
          if(isFALSE(on_ref_genome) && isTRUE(eff_on_minor_alleles)){
            stop("on_ref_genome must equal TRUE to use eff_on_minor_alleles")
          }
          if(isTRUE(on_ref_genome)){
            #Note this check will also work for when A0 is present
            switch_req <- get_genome_build(
              sumstats = copy(sumstats_dt),
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
                # Special case!! A0/A1 -> ref/alt so A1 flips meaning,
                # A0 is A* in mapping
                # but usually A1/A2 -> ref/alt so if A* found,
                # swap A1 to A2 and make A* -> A1
                new_headers <- colnames(sumstats_dt)
                if ("A*" %in% new_headers) {
                  # if A1 and A2 also present need to rename A2
                  if ("A1" %in% new_headers && "A2" %in% new_headers) {
                    data.table::setnames(sumstats_dt, "A2", "A2_from_input")
                  }
                  # if A1 present change to A2, doesn't have to be, 
                  # can be imputted
                  data.table::setnames(sumstats_dt, "A1", "A2", 
                                       skip_absent = TRUE)
                  data.table::setnames(sumstats_dt, "A*", "A1")
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
