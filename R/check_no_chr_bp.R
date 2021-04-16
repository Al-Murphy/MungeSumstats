#' Ensure that CHR and BP are missing if SNP is present, can find them
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS (GRCh37 or GRCh38)
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table setcolorder
#' @importFrom BSgenome snpsById
#' @importFrom data.table setorder
check_no_chr_bp <- function(sumstats_file, path, ref_genome){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = NULL
  # If SNP present but no CHR/BP then need to find them
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if(sum(c("CHR","BP") %in% col_headers)<=1 & sum("SNP" %in% col_headers)==1){
    SNP_LOC_DATA <- load_snp_loc_data(ref_genome,
                                        "Chromosome or Base Pair Position")
    sumstats_file <- data.table::fread(path)
    #if dataset has one of CHR or BP remove it and take from re dataset
    if(sum(c("CHR","BP") %in% col_headers)==1){
      colsToDelete <- c("CHR","BP")[c("CHR","BP") %in% col_headers]
      sumstats_file[,(colsToDelete):=NULL]
    }
    gr_rsids <- BSgenome::snpsById(SNP_LOC_DATA, ids = sumstats_file$SNP,
                                   ifnotfound="drop")#remove SNPs not found
    rsids <- data.table::setDT(data.frame(gr_rsids))
    data.table::setnames(rsids,"RefSNP_id","SNP")
    data.table::setorder(rsids,SNP)
    # join on CHR BP to sumstats
    data.table::setkey(sumstats_file,SNP)
    data.table::setkey(rsids,SNP)
    sumstats_file[rsids,CHR:=i.seqnames]
    sumstats_file[rsids,BP:=i.pos]
    #remove rows where CHR/BP couldn't be found
    sumstats_file <- sumstats_file[complete.cases(sumstats_file),]
    #move SNP, CHR, BP to start
    other_cols <-
      names(sumstats_file)[!names(sumstats_file) %in% c("SNP","CHR","BP")]
    data.table::setcolorder(sumstats_file, c("SNP","CHR","BP", other_cols))
    #write new data
    data.table::fwrite(sumstats_file,file=path,sep="\t")
    sumstats_file <- readLines(path)

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
