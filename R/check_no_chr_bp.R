#' Ensure that CHR and BP are missing if SNP is present, can find them
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return The modified sumstats_file
check_no_chr_bp <- function(sumstats_file){
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if(sum(c("CHR","BP") %in% col_headers)==0 & sum("SNP" %in% col_headers)==1){
    #sumstats = fread(path)
    #SNP_LOC_DATA = load_snp_loc_data()
    #SNP_LOC_DATA_2 = SNP_LOC_DATA[SNP_LOC_DATA$Build=="GRCh37",1:3]
    #sumstats2 = merge(sumstats,SNP_LOC_DATA_2,by="SNP")
    #sumstats3 = data.frame(sumstats2)[,c("SNP","CHR","BP",setdiff(colnames(sumstats2),c("SNP","CHR","BP")))]
    #fwrite(sumstats3,file=path,sep="\t"); sumstats_file <- readLines(path)
    stop("I've blocked this function because I've not tested it since replacing SNP_LOC_DATA. You should test it works manually. Let me know if you test it!")
  }#TO DO: add elseif version for either BP of CHR missing
  else(
    return(sumstats_file)
  )
}