#' Genome build liftover
#'
#' Transfer genomic coordinates from one genome build to another.
#'
#' @source \href{https://doi.org/doi:10.18129/B9.bioc.liftOver}{liftOver}
#' @source \href{https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/}{
#' UCSC chain files}
#'
#' @param sumstats_dt data table obj of the summary statistics
#'  file for the GWAS.
#' @param chrom_col Name of the chromosome column in 
#' \code{sumstats_dt} (e.g. "CHR").
#' @param start_col Name of the starting genomic position
#'  column in \code{sumstats_dt} (e.g. "POS","start").
#' @param end_col Name of the ending genomic position
#'  column in \code{sumstats_dt} (e.g. "POS","end"). 
#'  Can be the same as \code{start_col} when \code{sumstats_dt} 
#'  only contains SNPs that span 1 base pair (bp) each.
#' @param as_granges Return results as \link[GenomicRanges]{GRanges} 
#' instead of a \link[data.table]{data.table} (default: \code{FALSE}).
#' @param style Style to return \link[GenomicRanges]{GRanges} object in
#' (e.g.  "NCBI" = 4; "UCSC" = "chr4";) (default: \code{"NCBI"}).
#' @param verbose Print messages.
#' @inheritParams format_sumstats
#' 
#' @returns Lifted summary stats in \code{data.table} 
#' or \link[GenomicRanges]{GRanges} format.
#' 
#' @export
#' @importFrom rtracklayer liftOver width strand end
#' @importFrom GenomeInfoDb seqnames mapGenomeBuilds
#' @importFrom GenomicRanges mcols
#' @importFrom data.table as.data.table setnames :=
#' @examples 
#' sumstats_dt <- MungeSumstats::formatted_example()
#'
#' sumstats_dt_hg38 <- liftover(sumstats_dt=sumstats_dt, 
#'                              ref_genome = "hg19",
#'                              convert_ref_genome="hg38")
liftover <- function(sumstats_dt, 
                     convert_ref_genome, 
                     ref_genome, 
                     imputation_ind = TRUE,
                     chrom_col = "CHR",
                     start_col = "BP",
                     end_col = start_col, 
                     as_granges = FALSE,
                     style = "NCBI",
                     verbose = TRUE) {
    
    IMPUTATION_gen_build <- width <- strand <- end <- seqnames <- NULL;
    
    #### Map genome build synonyms ####
    query_ucsc <- if(!is.null(ref_genome)){
        GenomeInfoDb::mapGenomeBuilds(genome = ref_genome)$ucscID[1]
    } else {ref_genome}
    target_ucsc <- if(!is.null(convert_ref_genome)){
        GenomeInfoDb::mapGenomeBuilds(genome = convert_ref_genome)$ucscID[1]
    } else {convert_ref_genome}
     
    #### Check if one or more of the genomes couldn't be mapped ####
    null_builds <- c("query_genome", "target_genome")[
        c(is.null(query_ucsc), is.null(target_ucsc))
    ] 
    if(length(null_builds)>0){
        msg <- paste0("Could not recognize genome build of:\n",
                      paste(" -",null_builds,collapse = "\n"),
                      "\nThese will be inferred from the data.")
        message(msg)
    } 
    
    #### Check if liftover is necessary ####
    ## i.e. the desired genome build isn't the current one
    if ((!is.null(query_ucsc) & !is.null(target_ucsc)) &&
        (query_ucsc != target_ucsc)) {
        msg <- paste0(
            "Performing data liftover from ", query_ucsc, " to ",
            target_ucsc, "."
        )
        message(msg)

        #### Check that liftover is available ####
        ## If one or more builds are NULL, this won't be evaluated bc
        ## the builds will be inferred instead.
        if(query_ucsc=="hg38" && target_ucsc=="hg19") {
            build_conversion <- "hg38ToHg19"
        } else if (query_ucsc=="hg19" && target_ucsc=="hg38"){
            build_conversion <- "hg19ToHg38"
        } else {
            stop("Can only perform liftover between hg19 <---> hg38")
        } 
        #### Convert to GRanges #### 
        gr <- to_granges(
            sumstats_dt = sumstats_dt,
            style = "UCSC", 
            seqnames.field = chrom_col,
            start.field = start_col,
            end.field = end_col
        )
        #### Specify chain file ####
        chain <- get_chain_file(
            build_conversion = build_conversion,
            verbose = verbose
        )
        #### Liftover ####
        gr_lifted <- unlist(rtracklayer::liftOver(
            x = gr,
            chain = chain
        ))
        #### Chrom style ####
        gr_lifted <- granges_style(
            gr = gr_lifted,
            style = style
        )
        #### Return format ####
        if (as_granges) {
            if(imputation_ind){
                GenomicRanges::mcols(gr_lifted)[["IMPUTATION_gen_build"]] <-
                    TRUE
            }
            return(gr_lifted)
        } else {
            sumstats_dt <- data.table::as.data.table(gr_lifted)
            #### rename columns back to original #### 
            void_cols <- c("width","strand")
            void_cols <- void_cols[void_cols %in% names(sumstats_dt)] 
            if(length(void_cols)>0) sumstats_dt[,(void_cols):=NULL]
            data.table::setnames(sumstats_dt,"seqnames","CHR")
            #### Remove end_col if it was the same as start_col ####
            if (start_col == end_col) {
                sumstats_dt[, end := NULL]
                data.table::setnames(sumstats_dt, "start", start_col)
            } else {
                data.table::setnames(sumstats_dt, c("start","end"),
                                     c(start_col, end_col) 
                )
            }
            #### lastly rearrange the order again ####
            sumstats_dt <- check_col_order(
                sumstats_dt = sumstats_dt,
                path = NULL
            )$sumstats_dt
            if (imputation_ind) {
                sumstats_dt[, IMPUTATION_gen_build := TRUE]
            }
        } 
    } 
    return(sumstats_dt)
}
