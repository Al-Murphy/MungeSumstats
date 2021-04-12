#' Ensure that SNP is present if not can find it with CHR and BP
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table :=
#' @importFrom data.table setcolorder
#' @importFrom BSgenome snpsByOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
check_no_snp <- function(sumstats_file, path){
  SNP = CHR = i.RefSNP_id = NULL
  # If CHR and BP are present BUT not SNP then need to find the relevant SNP ids
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if(sum(c("CHR","BP") %in% col_headers)==2 & sum("SNP" %in% col_headers)==0){
    SNP_LOC_DATA <- load_snp_loc_data("SNP")
    sumstats_file <- data.table::fread(path)
    gr_snp <- GenomicRanges::makeGRangesFromDataFrame(sumstats_file,
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
    sumstats_file[,CHR:=tolower(as.character(CHR))]
    # join on SNP ID to sumstats
    data.table::setkeyv(sumstats_file,c("CHR","BP"))
    data.table::setkeyv(rsids,c("CHR","BP"))
    sumstats_file[rsids,SNP:=i.RefSNP_id]
    #move SNP to start
    other_cols <- names(sumstats_file)[names(sumstats_file)!="SNP"]
    data.table::setcolorder(sumstats_file, c("SNP", other_cols))
    #write new data
    data.table::fwrite(sumstats_file,file=path,sep="\t")
    sumstats_file <- readLines(path)

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
