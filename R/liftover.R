#' Genome build liftover
#'
#' Transfer genomic coordinates from one genome build to another.
#'
#' @source \href{https://doi.org/doi:10.18129/B9.bioc.liftOver}{liftOver}
#' @source \href{https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/}{
#' UCSC chain files}
#'
#' @inheritParams format_sumstats
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param verbose Print messages.
#' @return \code{list("sumstats_dt"=sumstats_dt)}
#' @keywords internal
#' @importFrom rtracklayer liftOver width strand end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom data.table data.table
liftover <- function(sumstats_dt, 
                     convert_ref_genome, 
                     ref_genome, 
                     imputation_ind,
                     verbose = TRUE) {
    IMPUTATION_gen_build <- NULL
    # check it's necessary i.e. the desired ref genome isn't the current one
    if (!is.null(convert_ref_genome) &&
        toupper(convert_ref_genome) != toupper(ref_genome)) {
        msg <- paste0(
            "Performing data liftover from ", ref_genome, " to ",
            convert_ref_genome, "."
        )
        message(msg)

        if (toupper(ref_genome) == "GRCH38") {#convert_ref_genome
            build_conversion <- "hg38ToHg19"
            ucsc_ref <- "hg38"
        } else {
            build_conversion <- "hg19ToHg38"
            ucsc_ref <- "hg19"
        }

        #### Convert to GRanges ####
        gr <- dt_to_granges(
            dat = sumstats_dt,
            style = "UCSC",
            chrom_col = "CHR",
            start_col = "BP",
            end_col = "BP"
        )
        #### Specify chain file ####
        chain <- get_chain_file(
            build_conversion = build_conversion,
            ucsc_ref = ucsc_ref,
            verbose = verbose
        )
        #### Liftover ####
        gr_lifted <- unlist(rtracklayer::liftOver(
            x = gr,
            chain = chain
        ))
        sumstats_dt <- as.data.table(gr_lifted)
        # rename columns back to org
        sumstats_dt[, width := NULL]
        sumstats_dt[, strand := NULL]
        sumstats_dt[, end := NULL]
        sumstats_dt[, seqnames := NULL]
        setnames(sumstats_dt, "start", "BP")
        # lastly rearrange the order again
        sumstats_dt <- check_col_order(
            sumstats_dt = sumstats_dt,
            path = NULL
        )$sumstats_dt
        if (imputation_ind) {
            sumstats_dt[, IMPUTATION_gen_build := TRUE]
        }
    }
    return(sumstats_dt)
}
