#' Ensure all SNPs are on the reference genome
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38). Default is GRCh37.
#' @param on_ref_genome Binary Should a check take place that all SNPs are on the reference genome by SNP ID. Default is TRUE
#' @param rsids datatable of snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' @return The modified sumstats_file
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table copy
check_on_ref_genome <- 
  function(sumstats_file, path, ref_genome, on_ref_genome, rsids){
  CHR = SNP = NULL
  # If CHR present and user specified chromosome to have SNPs removed
  col_headers <- names(sumstats_file)
  if("SNP" %in% col_headers && !isFALSE(on_ref_genome)){
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_file$SNP), ref_genome)
      #Save to parent environment so don't have to load again
      assign("rsids", rsids, envir = parent.frame())
    }
    num_bad_ids <- length(sumstats_file$SNP) - length(rsids$SNP)
    #check for SNPs not on ref genome
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs are not on the reference genome. ",
                    " These will be removed")
      message(msg)
      # join using SNP
      data.table::setkey(sumstats_file,SNP)
      #remove rows missing from the reference genome
      sumstats_file <- sumstats_file[rsids$SNP,]
    }
  }
  return(sumstats_file)
}
