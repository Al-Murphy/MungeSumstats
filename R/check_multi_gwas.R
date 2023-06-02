#' Ensure that only one model in GWAS sumstats or only one trait tested
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param analysis_trait If multiple traits were studied, name of the trait for
#' analysis from the GWAS. Default is NULL
#' @param mapping_file MungeSumstats has a pre-defined column-name mapping file
#' which should cover the most common column headers and their interpretations.
#' However, if a column header that is in youf file is missing of the mapping we
#' give is incorrect you can supply your own mapping file. Must be a 2 column
#' dataframe with column names "Uncorrected" and "Corrected". See
#' data(sumstatsColHeaders) for default mapping and necessary format.
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object
#' @keywords internal
check_multi_gwas <- function(sumstats_dt, path, analysis_trait, 
                             ignore_multi_trait,mapping_file) {
    if(isFALSE(ignore_multi_trait)){
      message("Checking for multi-GWAS.")
      col_headers <- names(sumstats_dt)
      # load synonym mapping - internal data no loading
      # Check for multiple p values and other columns,
      # if found choose 1 to continue
      # Use P value to find - first load all p-value synonyms
      p_syn <- sumstatsColHeaders[sumstatsColHeaders$Corrected == "P", 1]
      # So P_value columns found have already been corrected to "P" so if you find
      # any in p_syn in colheaders then we found multiple
      # NOTE THIS ASSUMES AN '_' or '.' SEPARATES P VALUE SYN FROM TRAIT/MODEL
      # NAME AND COLUMN NAME STARTS WITH P VAULE SYM
      # Sort p_syn so longest ones checked first, once 1 found stop
      # so don't get wrong subset with shorter ones
      p_syn <- p_syn[order(nchar(p_syn), p_syn, decreasing = TRUE)]
      for (pattern_i in p_syn) {
          p_val_multi <- grepl(paste0(
              "^", pattern_i, "_", "|",
              "^", pattern_i, "[.]"
          ), col_headers)
          if (any(p_val_multi)) {
              break
          }
      }
      # So note we could have 0 (do nothing), 1 (just update names),
      # 2+ (remove all but one) traits
      if (any(p_val_multi)) { # multiple traits
          # first get trait name by removing specific p_syn
          traits <- gsub(paste0(
              "^", pattern_i, "_", "|",
              "^", pattern_i, "[.]"
          ), "", col_headers[p_val_multi])
          if (length(traits) > 1) {
              if (!is.null(analysis_trait)) {
                  analysis_trait <- as.character(analysis_trait)
              } # ensure character
              msg2 <- paste0(
                  "WARNING: Multiple traits found in sumstats file only one",
                  " of which can be analysed: \n",
                  paste(traits, collapse = ", ")
              )
              message(msg2)
              stp_msg <- paste0(
                  "Inputted analysis_trait not one of these traits. Pleas",
                  "e input one of the traits above.\n  Or set",
                  "`ignore_multi_trait=TRUE` to ignore these multi-trait ",
                  "p-values and rerun `format_sumstats()`."
              )
              if (is.null(analysis_trait) ||
                  !toupper(analysis_trait) %in% toupper(traits)) {
                  stop(stp_msg)
              }
              # get right case for choice
              chosen_trait <- traits[toupper(traits) %in% toupper(analysis_trait)]
          } else { # just one trait
              msg <- paste0(
                "WARNING: A single multi-trait was found in sumstats file: \n",
                paste(traits, collapse = ", "),"\n",
                "If you don't want to use this value for P, set ",
                "`ignore_multi_trait=TRUE` and rerun `format_sumstats()`"
              )
              message(msg)
              chosen_trait <- traits
          }
          # just 1 trait from start or 1 chosen by user,
          # remove the trait from col name
          # then get the correct synonyms
          # NOTE ASSUMES TRAIT NAME STARTS WITH '_' or '.' &
          # IS AT THE END OF COLUMN
          chnge_header_names <- col_headers[grepl(
              paste0(
                  "_", chosen_trait, "$", "|",
                  "[.]", chosen_trait, "$"
              ),
              col_headers
          )]
  
          new_names <- gsub(paste0(
              "_", chosen_trait, "$", "|",
              ".", chosen_trait, "$"
          ), "", chnge_header_names)
          #if names present already, rename
          for(new_i in new_names){
              if(new_i %in% colnames(sumstats_dt)){
                  data.table::setnames(sumstats_dt, new_i, paste0(new_i,"_input"))
              }
          }
          data.table::setnames(sumstats_dt, chnge_header_names, new_names)
          # RE-standardise headers for all OS
          sumstats_return <-
              standardise_sumstats_column_headers_crossplatform(
                  sumstats_dt = sumstats_dt,
                  mapping_file = mapping_file
              )
          return(sumstats_return)
      } else {
          return(list("sumstats_dt" = sumstats_dt))
      }
    } else {
      return(list("sumstats_dt" = sumstats_dt))
    }
}
