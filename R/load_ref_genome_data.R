#' Load the reference genome data for SNPs of interest
#'
#' @param snps Character vector SNPs by rs_id from sumstats file of interest.
#' @param ref_genome Name of the reference genome used for the GWAS
#'  (GRCh37 or GRCh38)
#' @param msg Optional name of the column missing from the dataset in question.
#'  Default is NULL
#' @return data table of snpsById, filtered to SNPs of interest.
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table setorder
#' @importFrom BSgenome snpsById
load_ref_genome_data <- function(snps, ref_genome, msg = NULL) {
    message("Loading reference genome data.")
    SNP <- NULL
    SNP_LOC_DATA <- load_snp_loc_data(ref_genome, msg = msg)
    # Get correct ref genome
    if (toupper(ref_genome) == "GRCH37") {
        genome <-
            BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
    } else { # =="GRCH38"
        genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    }

    #snps can contain odd characters, need to sort these
    #so the snp should be all numeric of should start with rs and then be 
    #numeric. Find any that aren't and exclude these. The rs id can be inferred]
    #in a later function
    snps <-
        snps[!unlist(lapply(suppressWarnings(as.numeric(substring(snps,
                                                                  3,length(snps)
                                                                  ))),is.na))]
    gr_rsids <- BSgenome::snpsById(
        x = SNP_LOC_DATA,
        id = snps,
        genome = genome,
        ifnotfound = "drop"
    )
    rsids <- data.table::setDT(data.frame(gr_rsids))
    data.table::setnames(rsids, "RefSNP_id", "SNP")
    data.table::setorder(rsids, SNP)
    data.table::setkey(rsids, SNP)
    return(rsids)
}
