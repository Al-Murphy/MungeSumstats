#' Ensure that CHR and BP are missing if SNP is present, can find them
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
#' @importFrom BSgenome snpsById
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_no_chr_bp <- function(sumstats_dt, path, ref_genome,rsids){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = NULL
  # If SNP present but no CHR/BP then need to find them
  col_headers <- names(sumstats_dt)
  if(sum(c("CHR","BP") %in% col_headers)<=1 & sum("SNP" %in% col_headers)==1){
    #if dataset has one of CHR or BP remove it and take from re dataset
    if(sum(c("CHR","BP") %in% col_headers)==1){
      colsToDelete <- c("CHR","BP")[c("CHR","BP") %in% col_headers]
      sumstats_dt[,(colsToDelete):=NULL]
    }
    #check if rsids loaded if not do so
    if(is.null(rsids)){
      rsids <- load_ref_genome_data(copy(sumstats_dt$SNP), ref_genome,
                                      "Chromosome or Base Pair Position")
    }
    else{
      print_msg <- paste0("There is no Chromosome or Base Pair Position column",
                            " found within the data. It must be inferred from ",
                            " other column information.")
      message(print_msg)
    }
    # join on CHR BP to sumstats
    data.table::setkey(sumstats_dt,SNP)
    sumstats_dt[rsids,CHR:=i.seqnames]
    sumstats_dt[rsids,BP:=i.pos]
    #remove rows where CHR/BP couldn't be found
    sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt),]
    #move SNP, CHR, BP to start
    other_cols <-
      names(sumstats_dt)[!names(sumstats_dt) %in% c("SNP","CHR","BP")]
    data.table::setcolorder(sumstats_dt, c("SNP","CHR","BP", other_cols))

    return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt,"rsids"=rsids))
  }
}
