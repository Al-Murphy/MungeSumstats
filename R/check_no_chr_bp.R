#' Ensure that CHR and BP are missing if SNP is present, can find them
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param rsids datatable of snpsById, filtered to SNPs of interest if loaded already. Or else NULL
#' @return The modified sumstats_file
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table setcolorder
#' @importFrom BSgenome snpsById
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_no_chr_bp <- function(sumstats_file, path, ref_genome,rsids){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = NULL
  # If SNP present but no CHR/BP then need to find them
  col_headers <- names(sumstats_file)
  if(sum(c("CHR","BP") %in% col_headers)<=1 & sum("SNP" %in% col_headers)==1){
    #if dataset has one of CHR or BP remove it and take from re dataset
    if(sum(c("CHR","BP") %in% col_headers)==1){
      colsToDelete <- c("CHR","BP")[c("CHR","BP") %in% col_headers]
      sumstats_file[,(colsToDelete):=NULL]
    }
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_file$SNP), ref_genome,
                                      "Chromosome or Base Pair Position")
      #Save to parent environment so don't have to load again
      assign("rsids", rsids, envir = parent.frame())
    }
    else{
      print_msg <- paste0("There is no Chromosome or Base Pair Position column",
                            " found within the data. It must be inferred from ",
                            " other column information.")
      message(print_msg)
    }
    # join on CHR BP to sumstats
    data.table::setkey(sumstats_file,SNP)
    sumstats_file[rsids,CHR:=i.seqnames]
    sumstats_file[rsids,BP:=i.pos]
    #remove rows where CHR/BP couldn't be found
    sumstats_file <- sumstats_file[complete.cases(sumstats_file),]
    #move SNP, CHR, BP to start
    other_cols <-
      names(sumstats_file)[!names(sumstats_file) %in% c("SNP","CHR","BP")]
    data.table::setcolorder(sumstats_file, c("SNP","CHR","BP", other_cols))

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
