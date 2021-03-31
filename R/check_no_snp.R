#' Ensure that SNP is present if not can find it with CHR and BP
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
check_no_snp <- function(sumstats_file, path){
  # If CHR and BP are present BUT not SNP then need to find the relevant SNP ids
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]) 
  col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  
  if(sum(c("CHR","BP") %in% col_headers)==2 & sum("SNP" %in% col_headers)==0){
    print(paste0("There is no SNP column found within the data. ",
                  "It must be inferred from CHR and BP information."))
    msg<-paste0("Which genome build is the data from? ",
                  "1 for GRCh37, 2 for GRCh38: ")
    genomebuild <- as.numeric(readline(msg))
    if(!genomebuild %in% c(1,2))
      stop(paste0("Genome build must be entered as either ",
                    "1 (for GRCh37) or 2 (for GRCh38)"))
    
    SNP_LOC_DATA = load_snp_loc_data()
    if(genomebuild==1){
      genomebuild="GRCh37"
    }
    else{
      genomebuild="GRCh38"
    }
    snpLocDat = SNP_LOC_DATA[SNP_LOC_DATA$Build==genomebuild,][,-4]
    #Use data table for speed
    sumstats = fread(path)
    sumstats$CHR = as.factor(sumstats$CHR)
    if(length(grep("chr",sumstats$CHR[1]))!=0)
      sumstats$CHR = gsub("chr","",sumstats$CHR)
    sumstats2 <- merge(sumstats,snpLocDat,by=c("CHR","BP"))
    # Remove any rows where P is NaN
    sumstats2 <- sumstats2[!is.nan(sumstats2$P),]
    fwrite(sumstats2,file=path,sep="\t")
    sumstats_file <- readLines(path)
    
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}