#' Ensure that SNP is present if not can find it with CHR and BP
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @param verbose should messages be printed. Default it TRUE.
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table :=
#' @importFrom data.table setcolorder
#' @importFrom data.table copy
#' @importFrom BSgenome snpsByOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
check_no_snp <- function(sumstats_dt, path, ref_genome, verbose=TRUE){
  SNP = CHR = i.RefSNP_id = NULL
  # If CHR and BP are present BUT not SNP then need to find the relevant SNP ids
  col_headers <- names(sumstats_dt)
  if(sum(c("CHR","BP") %in% col_headers)==2 & sum("SNP" %in% col_headers)==0){
    msg <- "SNP"
    if(isFALSE(verbose))
      msg <- NULL
    SNP_LOC_DATA <- load_snp_loc_data(ref_genome,msg)
    #if chromosome col has chr prefix remove it
    sumstats_dt[,CHR:=gsub("chr","",CHR)]
    gr_snp <- GenomicRanges::makeGRangesFromDataFrame(copy(sumstats_dt),
                                                      keep.extra.columns = TRUE,
                                                      seqnames.field = "CHR",
                                                      start.field = "BP",
                                                      end.field = "BP")
    gr_rsids <-
      BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = gr_snp)
    rsids <- data.table::setDT(data.frame(gr_rsids))
    data.table::setnames(rsids,"seqnames","CHR")
    data.table::setnames(rsids,"pos","BP")
    #in case there is CHR8 and chr8
    rsids[,CHR:=tolower(as.character(CHR))]
    sumstats_dt[,CHR:=tolower(as.character(CHR))]
    # join on SNP ID to sumstats
    data.table::setkeyv(sumstats_dt,c("CHR","BP"))
    data.table::setkeyv(rsids,c("CHR","BP"))
    sumstats_dt[rsids,SNP:=i.RefSNP_id]
    #remove rows where SNP couldn't be found
    sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt),]
    #move SNP to start
    other_cols <- names(sumstats_dt)[names(sumstats_dt)!="SNP"]
    data.table::setcolorder(sumstats_dt, c("SNP", other_cols))

    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
