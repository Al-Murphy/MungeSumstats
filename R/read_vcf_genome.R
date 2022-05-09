#' Read VCF genome
#' 
#' Get the genome build of a remote or local VCF file.
#' @param header Header extracted by \link[VariantAnnotation]{scanVcfHeader}.
#' @param validate Validate genome name using 
#' \link[GenomeInfoDb]{mapGenomeBuilds}.
#' @param default_genome When no genome can be extracted,
#'  default to this genome build. 
#' @param verbose Print messages. 
#' 
#' @keywords internal
#' @importFrom GenomeInfoDb mapGenomeBuilds
#' @importFrom VariantAnnotation genome
#' @returns genome 
read_vcf_genome <- function(header=NULL,
                            validate=FALSE,
                            default_genome="HG19/GRCh37",
                            verbose=TRUE){
    
    msg <- paste("No known genome build matches.",
                 "Defaulting to:",default_genome)
    genome <- tryCatch({
        #### Try 1 ###
        genome <- as.character(VariantAnnotation::genome(header)[[1]])
        #### Try 2 ###
        if(length(genome)==0){
            genome <- as.character(header@header$contig$assembly[[1]])
        } 
        if(isTRUE(validate)){
            genome_matches <- GenomeInfoDb::mapGenomeBuilds(genome=genome)
            if(nrow(genome_matches)==0) {
                messager(msg, v=verbose)
                genome <- default_genome
            }
        } 
        genome
    }, error = function(e){messager(msg, v=verbose); default_genome})
    return(genome)
}
