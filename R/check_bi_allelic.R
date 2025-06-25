#' Remove non-biallelic SNPs
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics data table object
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest
#'   if loaded already. Or else \code{NULL}.
#'   \item \code{log_files}: log file list
#' }
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table copy
#' @importFrom Biostrings IUPAC_CODE_MAP
check_bi_allelic <-
    function(sumstats_dt, path, ref_genome, bi_allelic_filter, rsids,
             log_folder_ind, check_save_out, tabix_index, nThread, log_files,
             dbSNP,dbSNP_tarball) {
        CHR <- alleles_as_ambig <- SNP <- NULL
        # If SNP present and user specified to remove
        col_headers <- names(sumstats_dt)
        if ("SNP" %in% col_headers && !isFALSE(bi_allelic_filter)) {
            message("Checking for bi-allelic SNPs.")
            if (is.null(rsids)) {
                rsids <- load_ref_genome_data(
                    data.table::copy(sumstats_dt$SNP),
                    ref_genome,
                    dbSNP = dbSNP,
                    dbSNP_tarball = dbSNP_tarball
                )
            }
            # get chars for SNPs not bi/tri allelic
            # or strand ambig from IUPAC_CODE_MAP
            nonambig_IUPAC_CODE_MAP <-
                names(Biostrings::IUPAC_CODE_MAP[nchar(
                    Biostrings::IUPAC_CODE_MAP
                ) < 3])
            # ensure rsids is up-to-date with filtered sumstats_dt
            rsids <- rsids[unique(sumstats_dt$SNP), , nomatch = NULL]
            data.table::setkey(rsids, SNP)
            num_bad_ids <- nrow(
                rsids[!alleles_as_ambig %in% nonambig_IUPAC_CODE_MAP]
            )
            # check for SNPs not on ref genome
            if (num_bad_ids > 0) {
                msg <- paste0(
                    formatC(num_bad_ids, big.mark = ","),
                    " SNPs are non-biallelic.",
                    " These will be removed."
                )
                message(msg)
                # join using SNP
                data.table::setkey(sumstats_dt, SNP)
                keep_snps <- rsids[
                    alleles_as_ambig %in% nonambig_IUPAC_CODE_MAP
                ]$SNP
                # remove rows not bi-allelic
                # If user wants log, save it to there
                if (log_folder_ind) {
                    name <- "snp_bi_allelic"
                    name <- get_unique_name_log_file(
                        name = name,
                        log_files = log_files
                    )
                    write_sumstats(
                        sumstats_dt = sumstats_dt[!keep_snps, ],
                        save_path =
                            paste0(
                                check_save_out$log_folder,
                                "/", name,
                                check_save_out$extension
                            ),
                        sep = check_save_out$sep,
                        #don't tab indx as could be miss values & cause err
                        #tabix_index = tabix_index,
                        nThread = nThread
                    )
                    log_files[[name]] <-
                        paste0(
                            check_save_out$log_folder, "/", name,
                            check_save_out$extension
                        )
                }
                sumstats_dt <- sumstats_dt[keep_snps, ]
            }
            return(list(
                "sumstats_dt" = sumstats_dt, "rsids" = rsids,
                "log_files" = log_files
            ))
        } else {
            return(list(
                "sumstats_dt" = sumstats_dt, "rsids" = rsids,
                "log_files" = log_files
            ))
        }
    }
