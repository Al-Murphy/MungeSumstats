#' Loads the SNP locations and alleles for Homo sapiens extracted from
#' NCBI dbSNP Build 144. Reference genome version is dependent on user input.
#'
#' @param ref_genome name of the reference genome used for the GWAS 
#' (GRCh37 or GRCh38)
#' @param dbSNP version of dbSNP to be used (144 or 155)
#' @param msg Optional name of the column missing from the dataset in question
#' @return SNP_LOC_DATA SNP positions and alleles for Homo sapiens extracted
#' from NCBI dbSNP Build 144
#'
#' @examples
#' SNP_LOC_DATA <- load_snp_loc_data("GRCH37",dbSNP=144)
#' @export
load_snp_loc_data <- function(ref_genome,dbSNP=c(144,155),
                              msg = NULL) {
    ref_genome <- toupper(ref_genome)
    message("Loading SNPlocs data.")
    if (!is.null(msg)) {
        print_msg <- paste0(
            "There is no ", msg, " column found within the data. ",
            "It must be inferred from other column information."
        )
        message(print_msg)
    }
    if (ref_genome == "GRCH37") {
        if(dbSNP==144){
          snp_loc_data <-
              SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
        }else{ #155
          stop("dbSNP 155 for GRCh37 not implemented......")
          snp_loc_data <-
            SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37
        }
    } else { # =="GRCH38"
        if(dbSNP==144){
          snp_loc_data <-
              SNPlocs.Hsapiens.dbSNP144.GRCh38::SNPlocs.Hsapiens.dbSNP144.GRCh38
        } else{ #155
          snp_loc_data <-
            SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
        }  
    }
    return(snp_loc_data)
}
