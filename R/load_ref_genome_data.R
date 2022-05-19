#' Load the reference genome data for SNPs of interest
#'
#' @param snps Character vector SNPs by rs_id from sumstats file of interest.
#' @param ref_genome Name of the reference genome used for the GWAS
#'  (GRCh37 or GRCh38)
#' @param msg Optional name of the column missing from the dataset in question.
#'  Default is NULL
#' @return data table of snpsById, filtered to SNPs of interest.
#' @importFrom data.table setDT setkey := setnames setorder
#' @importFrom BSgenome snpsById
#' @importFrom stringr str_sub
#' @source 
#' \code{
#' sumstats_dt <- formatted_example()
#' rsids <- MungeSumstats:::load_ref_genome_data(snps = sumstats_dt$SNP,
#'                                               ref_genome = "GRCH37")
#' }
load_ref_genome_data <- function(snps, 
                                 ref_genome, 
                                 msg = NULL) {
    
    SNP <- NULL
    SNP_LOC_DATA <- load_snp_loc_data(ref_genome, msg = msg)
    # Get correct ref genome
    message("Loading reference genome data.")
    if (toupper(ref_genome) == "GRCH37") {
        genome <-
            BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
    } else { # =="GRCH38"
        genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    }

    #snps can contain odd characters, need to sort these
    #so the snp should be all numeric or should start with rs and then be 
    #numeric. Find any that aren't and exclude these. The rs id can be inferred
    #in a later function
    messager("Preprocessing RSIDs.") 
    snp_check <- suppressMessages(as.numeric(
        stringr::str_sub(string = snps, start = 3, end = -1L) 
    ))
    snps <- snps[!is.na(snp_check)]
    #need to make snps unique for cases where non bi-allelic or indels
    snps <- unique(snps)
    
    messager("Validating RSIDs of",formatC(length(snps),big.mark = ","),
             "SNPs using BSgenome::snpsById...")
    tm <- system.time({
        gr_rsids <- BSgenome::snpsById(
            x = SNP_LOC_DATA,
            id = snps,
            genome = genome,
            ifnotfound = "drop"
        )
    }) 
    messager("BSgenome::snpsById done in",round(tm[[3]]),"seconds.")
     
    rsids <- data.table::setDT(data.frame(gr_rsids))
    data.table::setnames(rsids, "RefSNP_id", "SNP")
    data.table::setorder(rsids, SNP)
    data.table::setkey(rsids, SNP)
    return(rsids)
}
