#' Remove SNPs with strand-ambiguous alleles
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
#' @importFrom data.table setkey
#' @importFrom data.table :=
check_strand_ambiguous <- function(sumstats_dt, 
                                   path, 
                                   ref_genome, 
                                   strand_ambig_filter,
                                   log_folder_ind,
                                   check_save_out, 
                                   tabix_index, 
                                   nThread, 
                                   log_files) {
        CHR <- alleles_as_ambig <- SNP <- A1 <- A2 <- NULL
        # If SNP present and user specified to remove
        col_headers <- names(sumstats_dt)
        if ("SNP" %in% col_headers && !isFALSE(strand_ambig_filter)) {
            message("Checking for strand ambiguous SNPs.")
            A_T_ambig <- sumstats_dt[A1 == "A" & 
                                         A2 == "T" | 
                                         A1 == "T" & 
                                         A2 == "A", ]$SNP
            C_G_ambig <- sumstats_dt[A1 == "C" &
                                         A2 == "G" | 
                                         A1 == "G" & 
                                         A2 == "C", ]$SNP
            num_bad_ids <- length(A_T_ambig) + length(C_G_ambig)
            # check for SNPs not on ref genome
            if (num_bad_ids > 0) {
                msg <- paste0(
                    formatC(num_bad_ids, big.mark = ","),
                    " SNPs are strand-ambiguous alleles including",
                    " ", formatC(length(A_T_ambig), big.mark = ","),
                    " A/T and ", formatC(length(C_G_ambig), big.mark = ","),
                    " C/G ambiguous SNPs. These will be removed"
                )
                message(msg)
                rmv_snps <- c(A_T_ambig, C_G_ambig)
                # join using SNP
                data.table::setkey(sumstats_dt, SNP)
                # remove strand ambiguous SNPs
                # If user wants log, save it to there
                if (log_folder_ind) {
                    name <- "snp_strand_ambiguous"
                    name <- get_unique_name_log_file(name = name,
                                                     log_files = log_files)
                    write_sumstats(
                        sumstats_dt = sumstats_dt[rmv_snps, ],
                        save_path =
                            paste0(
                                check_save_out$log_folder,
                                "/", name,
                                check_save_out$extension
                            ),
                        sep = check_save_out$sep,
                        tabix_index = tabix_index,
                        nThread = nThread
                    )
                    log_files[[name]] <-
                        paste0(
                            check_save_out$log_folder, "/", name,
                            check_save_out$extension
                        )
                }
                sumstats_dt <- sumstats_dt[!rmv_snps, ]
            }
            return(list("sumstats_dt" = sumstats_dt,
                        "log_files" = log_files))
        } else {
            return(list("sumstats_dt" = sumstats_dt,
                        "log_files" = log_files))
        }
    }
