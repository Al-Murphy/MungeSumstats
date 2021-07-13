#' Infers the genome build of the summary statistics file (GRCh37 or GRCh38) 
#' from the data
#'
#' @param sumstats data table/data frame obj of the summary statistics file for
#' the GWAS or file path to summary statistics file
#' @param nThread Number of threads to use for parallel processes. 
#' @return ref_genome the genome build of the data
#'
#' @examples
#' #Pass path to Educational Attainment Okbay sumstat file to a temp directory
#' 
#' eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt", 
#'                                  package="MungeSumstats")
#' ref_genome <- get_genome_build(eduAttainOkbayPth)                                  
#'
#' @export
#' @importFrom data.table setDT
get_genome_build <- function(sumstats, nThread = 1){
  #if not a data.table, must be a path
  if(!is.data.frame(sumstats)){
    # Checking if the file exists should happen first
    if (!file.exists(sumstats))
      stop("Path to GWAS sumstats is not valid") 
    #read in data
    sumstats <- read_sumstats(path = sumstats, 
                                nThread = nThread)
  }else{
    #ensure data table obj
    sumstats <- data.table::setDT(sumstats)
  }
  #need SNP ID column (RS ID) CHR and BP (POS) to infer build - check these are 
  #present, considering all known names
  sumstats_return <-
    standardise_sumstats_column_headers_crossplatform(sumstats_dt = sumstats,
                                                      path =  NULL)
  sumstats <- sumstats_return$sumstats_dt
  
  
  err_msg <- 
    paste0("SNP ID column (RS ID) CHR and BP (POSITION) columns are needed to ",
           "infer the genome build. These could not be found in your dataset. ",
           "Please specify the genome build manually to run format_sumstats()")
  if(!all(c("SNP","CHR","BP") %in% colnames(sumstats)))
    stop(err_msg)
  #otherwise SNP, CHR, BP were all found and can infer
  snp_loc_data_37 <- load_ref_genome_data(snps = sumstats$SNP,
                                            ref_genome = "GRCH37")
  snp_loc_data_38 <- load_ref_genome_data(snps = sumstats$SNP,
                                          ref_genome = "GRCH38")
  #convert CHR filed in ref genomes to character not factor
  snp_loc_data_37[,seqnames:=as.character(seqnames)]
  snp_loc_data_38[,seqnames:=as.character(seqnames)]
  #convert CHR filed in data to character if not already
  sumstats[,CHR:=as.character(CHR)]
  #Now check which genome build has more matches to data
  num_37 <- 
    nrow(snp_loc_data_37[sumstats, , on = c(SNP="SNP", pos="BP",seqnames="CHR"),
                          nomatch=FALSE])
  num_38 <- 
    nrow(snp_loc_data_37[sumstats, , on = c(SNP="SNP", pos="BP",seqnames="CHR"),
                         nomatch=FALSE])
  if(num_37>num_38){
    ref_gen_num <- num_37
    ref_genome <- "GRCH37"
  }  
  else{
    ref_gen_num <- num_38
    ref_genome <- "GRCH38"
  }  
  #add a warning if low proportion of matches found
  msg <- paste0("WARNING: Less than 10% of your SNPs matched that of either ",
                "reference genome, this question the quality of your summary ",
                "statistics file.")
  if(ref_gen_num/nrow(sumstats)<0.1)
    message(msg)
  return(ref_genome)
}
