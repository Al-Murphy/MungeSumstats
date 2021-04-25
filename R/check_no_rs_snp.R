#' Ensure that SNP starts with rs
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table :=
#' @importFrom data.table setcolorder
#' @importFrom data.table copy
#' @importFrom data.table tstrsplit
#' @importFrom data.table rbindlist
#' @importFrom BSgenome snpsByOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
check_no_rs_snp <- function(sumstats_dt, path, ref_genome){
  # test case has X:102670736:T_TAT, chr1:60320992 and rs2326918 all in SNP col   
  SNP = CHR = CHR1 = BP1 = i.RefSNP_id = NULL
  # If SNP column doesn't start with rs
  col_headers <- names(sumstats_dt)
  if(sum("SNP" %in% col_headers)==1){
    miss_rs <- sumstats_dt[!grep("^rs",SNP),]
    #first case is chr:bp together - impute SNP for these
    miss_rs_chr_bp <- miss_rs[grep(":",SNP),]
    if(nrow(miss_rs)!=nrow(sumstats_dt) && nrow(miss_rs_chr_bp)>0){
      #remove snps missing rs
      sumstats_dt <- sumstats_dt[grep("^rs",SNP),]
      #check if any have more than 1 ":" remove these
      miss_rs_chr_bp <- miss_rs_chr_bp[!grep(".*:.*:.*",SNP)]
      msg <- paste0(nrow(miss_rs_chr_bp)," SNP IDs appear to be made up of ",
                      "chr:bp, these will be replaced by their SNP ID from the",
                      " reference genome")
      message(msg)
      SNP_LOC_DATA <- load_snp_loc_data(ref_genome,NULL)
      #split out chr:bp - check if chr or bp first by longer of two
      splits <- strsplit(miss_rs_chr_bp[1,SNP],split=':', fixed=TRUE)[[1]]
      if(nchar(splits[1])<nchar(splits[2]))
        format <- c("CHR1","BP1")
      else
        format <- c("BP1","CHR1")
      miss_rs_chr_bp[, (format) := data.table::tstrsplit(SNP,
                                                      split=":", fixed=TRUE)]
      #ensure integer col
      miss_rs_chr_bp[,BP1:=as.integer(BP1)]
      #now drop SNP
      miss_rs_chr_bp[,SNP := NULL]
      #if chromosome col has chr prefix remove it
      miss_rs_chr_bp[,CHR1:=gsub("chr","",CHR1)]
      gr_snp <- 
        GenomicRanges::makeGRangesFromDataFrame(copy(miss_rs_chr_bp),
                                                keep.extra.columns = TRUE,
                                                seqnames.field = "CHR1",
                                                start.field = "BP1",
                                                end.field = "BP1")
      gr_rsids <-
        BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = gr_snp)
      rsids <- data.table::setDT(data.frame(gr_rsids))
      data.table::setnames(rsids,"seqnames","CHR1")
      data.table::setnames(rsids,"pos","BP1")
      #in case there is CHR8 and chr8
      rsids[,CHR1:=tolower(as.character(CHR1))]
      miss_rs_chr_bp[,CHR1:=tolower(as.character(CHR1))]
      # join on SNP ID to sumstats
      data.table::setkeyv(miss_rs_chr_bp,c("CHR1","BP1"))
      data.table::setkeyv(rsids,c("CHR1","BP1"))
      miss_rs_chr_bp[rsids,SNP:=i.RefSNP_id]
      #remove rows where SNP couldn't be found
      miss_rs_chr_bp <- miss_rs_chr_bp[complete.cases(miss_rs_chr_bp),]
      #remove temp colunmns
      miss_rs_chr_bp[,(format):=NULL]
      #get columns in same order as rest of data table
      data.table::setcolorder(miss_rs_chr_bp, col_headers)
      #join with full dataset
      sumstats_dt <- data.table::rbindlist(list(sumstats_dt,miss_rs_chr_bp))
    }
    if(nrow(miss_rs)!=nrow(sumstats_dt) && nrow(miss_rs)!=0){
      if(nrow(miss_rs_chr_bp)==0){#don't filter twice if hit prev condition
        #remove snps missing rs
        sumstats_dt <- sumstats_dt[grep("^rs",SNP),]
      }  
      #if any weird SNP rows left that aren't chr:bp or rs id's remove them
      msg <- paste0(nrow(miss_rs) - nrow(miss_rs_chr_bp)," SNP IDs are not ",
                    "correctly formatted and will be removed")
      message(msg)
    }
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
