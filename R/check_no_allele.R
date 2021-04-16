#' Ensure that A1 & A2 are present, if not can find it with SNP and other allele
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
#' @importFrom data.table setorder
#' @importFrom BSgenome snpsById
#' @importFrom BSgenome.Hsapiens.NCBI.GRCh38 BSgenome.Hsapiens.NCBI.GRCh38
#' @importFrom BSgenome.Hsapiens.1000genomes.hs37d5 BSgenome.Hsapiens.1000genomes.hs37d5
check_no_allele <- function(sumstats_file, path, ref_genome){
  SNP = i.seqnames = CHR = BP = i.pos = LP = P = A1 = A2 = 
    i.ref_allele = i.alt_alleles = alt_alleles =NULL
  # If SNP present but no A1/A2 then need to find them
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if(sum(c("A1","A2") %in% col_headers)<=1 & sum("SNP" %in% col_headers)==1){
    #save(sumstats_file, path, ref_genome,file="~/Downloads/temp.rda")
    #print(aa)
    SNP_LOC_DATA <- load_snp_loc_data(ref_genome,"A1 or A2")
    sumstats_file <- data.table::fread(path)
    #Get correct ref genome
    if(toupper(ref_genome)=="GRCH37"){
      genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    }
    else{#=="GRCH38"
      genome <- BSgenome.Hsapiens.NCBI.GRCh38
    }
    gr_rsids <- snpsById(SNP_LOC_DATA, sumstats_file$SNP, 
                          genome=genome,ifnotfound="drop")
    #NOTE there is only one reference allele (A1) but
    #there can be multiple alternate alleles (A2) but most just have one
    #so if dataset has reference allele (A1) just remove it and 
    #get it from ref genome
    if(sum(c("A1") %in% col_headers)==1)
      sumstats_file[,A1:=NULL]
    
    rsids <- data.table::setDT(data.frame(gr_rsids))
    data.table::setnames(rsids,"RefSNP_id","SNP")
    data.table::setorder(rsids,SNP)
    # join based on SNP as key
    data.table::setkey(sumstats_file,SNP)
    data.table::setkey(rsids,SNP)
    
    #so if dataset has alternate alleles (A2), use it rather than ref genomes 
    #get it from ref genome
    if(sum(c("A2") %in% col_headers)==1){
      sumstats_file[rsids,A1:=i.ref_allele]
    }
    else{ #get both A1, A2 from ref genome - choose an A2 value where multiple
      sumstats_file[rsids,A1:=i.ref_allele]
      msg <- paste0("WARNING: Inferring the alternative allele (A2) from the ",
                    "reference genome. In some instances, there are more ",
                    "than one\nalternative allele. Arbitrarily, only the ",
                    "first will be kept. See column `alt_alleles` in your ",
                    "returned sumstats file\nfor all alternative alleles.")
      message(msg)
      sumstats_file[rsids,alt_alleles:=i.alt_alleles]
      #just take first A2 value arbitrarily
      sumstats_file[,A2:= as.character(lapply(alt_alleles, function(x) x[1]))]
      #collapse alt_alleles into character type sep by columns
      sumstats_file[,alt_alleles:= 
                      as.character(lapply(alt_alleles, 
                                          function(x) paste0(x,
                                                             collapse = ",")))]
    }
    #remove rows where A1/A2 couldn't be found
    sumstats_file <- sumstats_file[complete.cases(sumstats_file),]
    #move SNP, CHR, BP, A1 and A2 to start
    other_cols <-
      names(sumstats_file)[!names(sumstats_file) %in% 
                             c("SNP","CHR","BP", "A1", "A2")]
    data.table::setcolorder(sumstats_file, 
                              c("SNP","CHR","BP", "A1", "A2", other_cols))
    #write new data
    data.table::fwrite(sumstats_file,file=path,sep="\t")
    sumstats_file <- readLines(path)
    
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
