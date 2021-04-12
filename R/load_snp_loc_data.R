#' Loads the SNP locations and alleles for Homo sapiens extracted from
#' NCBI dbSNP Build 144. Reference genome version is dependent on user input.
#'
#' @param msg Optional name of the column missing from the dataset in question
#' @return SNP_LOC_DATA SNP positions and alleles for Homo sapiens extracted
#' from NCBI dbSNP Build 144
#'
#' @examples
#' #SNP_LOC_DATA <- load_snp_loc_data()
#'
#' @export
#' @importFrom SNPlocs.Hsapiens.dbSNP144.GRCh37 SNPlocs.Hsapiens.dbSNP144.GRCh37
#' @importFrom SNPlocs.Hsapiens.dbSNP144.GRCh38 SNPlocs.Hsapiens.dbSNP144.GRCh38
load_snp_loc_data <- function(msg=NULL){
  if(!is.null(msg)){
    message(paste0("There is no ",msg," column found within the data. ",
                    "It must be inferred from other column information."))
  }
  genomebuild <- 0
  while((!genomebuild %in% c(1,2))){
    msg<-paste0("Which genome build is the data from? ",
                "1 for GRCh37, 2 for GRCh38: ")
    genomebuild <- as.numeric(readline(msg))
    if(!genomebuild %in% c(1,2))
      message(paste0("Genome build must be entered as either ",
                      "1 (for GRCh37) or 2 (for GRCh38)"))
  }
  if(genomebuild==1){
    snp_loc_data <-
      SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
  }
  else{
    snp_loc_data <-
      SNPlocs.Hsapiens.dbSNP144.GRCh38::SNPlocs.Hsapiens.dbSNP144.GRCh38
  }
  return(snp_loc_data)
}
