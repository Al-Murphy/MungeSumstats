#' Ensure that A1 & A2 are present, if not can find it with SNP and other allele
#'
#' More care needs to be taken if one of A1/A2 is present, before imputing the
#' other allele flipping needs to be checked
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics data table object
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest 
#'   if loaded already. Or else NULL.
#'   \item \code{allele_flip_check}: does the dataset require allele flip check
#'   \item \code{log_files}: log file list
#'   \item \code{bi_allelic_filter}: should multi-allelic SNPs be filtered out
#' }
#' @keywords internal
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table setcolorder
#' @importFrom data.table setorder
#' @importFrom data.table copy
#' @importFrom Biostrings IUPAC_CODE_MAP
check_no_allele <- function(sumstats_dt, 
                            path, 
                            ref_genome, 
                            rsids,
                            imputation_ind, 
                            allele_flip_check, 
                            log_folder_ind,
                            check_save_out, 
                            tabix_index, 
                            nThread, 
                            log_files,
                            bi_allelic_filter,
                            dbSNP) {
    SNP <- i.seqnames <- CHR <- BP <- i.pos <- LP <- P <- A1 <- A2 <-
        i.ref_allele <- i.alt_alleles <- alt_alleles <- ss_A <-
        alleles_as_ambig <- IMPUTATION_A2 <- IMPUTATION_A1 <- NULL
    # If SNP present but no A1/A2 then need to find them
    col_headers <- names(sumstats_dt)
    if (sum(c("A1", "A2") %in% col_headers) <= 1 &
        sum("SNP" %in% col_headers) == 1) {
        # if A2 is missing need to remove multi-allelic SNPs 
        # as we don't know which is the one the author is thinking of
        msg_a2 <- paste0(
            "WARNING: No A2 column found in the data, multi-allelic ",
            "can't not be accurately chosen (as any\nof the choices ",
            "could be valid). bi_allelic_filter has been forced to ",
            "TRUE."
        )
        if (!"A2" %in% col_headers && !bi_allelic_filter) {
            bi_allelic_filter <- TRUE
            message(msg_a2)
        }
        # check if rsids loaded if not do so
        if (is.null(rsids)) {
            rsids <- load_ref_genome_data(
                data.table::copy(sumstats_dt$SNP),
                ref_genome = ref_genome,
                dbSNP = dbSNP,
                msg = "A1 or A2 allele information"
            )
        } else {
            print_msg <- paste0(
                "There is no A1 or A2 allele information column",
                " found within the data. It must be inferred from ",
                " other column information."
            )
            message(print_msg)
        }
        # ensure rsids is up-to-date with filtered sumstats_dt
        rsids <- rsids[unique(sumstats_dt$SNP), , nomatch = NULL]
        data.table::setkey(rsids, SNP)

        # join based on SNP as key
        data.table::setkey(sumstats_dt, SNP)
        # if one allele in dataset join other
        if (sum(c("A1", "A2") %in% col_headers) == 1) { # one allele missing
            # NOTE: This is the old flipping code but doesn't matter since the
            # flipping check will correct
            message("One of A1/A2 are missing, allele flipping will be tested")
            allele_flip_check <- TRUE
            msg_a <- paste0(
                "WARNING: One of A1/A2 are missing, bi_allelic_filter ",
                "has been forced to TRUE to test imputation."
            )
            if (!allele_flip_check) {
                allele_flip_check <- TRUE
                message(msg_a)
            }
            # First get col name in data for allele we have
            ssA <- c("A1", "A2")[c("A1", "A2") %in% col_headers]
            # join based on SNP as key
            data.table::setorder(sumstats_dt, SNP)
            data.table::setkey(sumstats_dt, SNP)
            # different rules for each
            if ("A1" == ssA) {
                message("Deriving A2 from reference genome")
                msg <- paste0(
                    "WARNING: Inferring the alternative allele (A2) from the",
                    " reference genome. In some instances, there are more ",
                    "than one\nalternative allele. Arbitrarily, only the ",
                    "first will be kept. See column `alt_alleles` in your ",
                    "returned sumstats file\nfor all alternative alleles."
                )
                message(msg)
                sumstats_dt[rsids, alt_alleles := i.alt_alleles]
                # just take first A2 value arbitrarily
                sumstats_dt[, A2 := as.character(
                    lapply(alt_alleles,function(x) x[1]))]
                # collapse alt_alleles into character type sep by columns
                sumstats_dt[, alt_alleles :=
                    as.character(lapply(
                        alt_alleles,
                        function(x) {
                            paste0(x,
                                collapse = ","
                            )
                        }
                    ))]
                # if user specifies add a column to notify of the imputation
                if (imputation_ind) {
                    sumstats_dt[, IMPUTATION_A2 := TRUE]
                }
            } else { # A2 == ssA
                message("Deriving A1 from reference genome")
                sumstats_dt[rsids, A1 := i.ref_allele]
                # if user specifies add a column to notify of the imputation
                if (imputation_ind) {
                    sumstats_dt[, IMPUTATION_A1 := TRUE]
                }
            }
        } else {
            # get both A1, A2 from ref genome - 
            # choose an A2 value where multiple
            message("Deriving both A1 and A2 from reference genome")
            sumstats_dt[rsids, A1 := i.ref_allele]
            msg <- paste0(
                "WARNING: Inferring the alternative allele (A2) from the ",
                "reference genome. In some instances, there are more ",
                "than one\nalternative allele. Arbitrarily, only the ",
                "first will be kept. See column `alt_alleles` in your ",
                "returned sumstats file\nfor all alternative alleles."
            )
            message(msg)
            sumstats_dt[rsids, alt_alleles := i.alt_alleles]
            # just take first A2 value arbitrarily
            sumstats_dt[, A2 := as.character(
                lapply(alt_alleles, function(x) x[1]))]
            # collapse alt_alleles into character type sep by columns
            sumstats_dt[, alt_alleles :=
                as.character(lapply(
                    alt_alleles,
                    function(x) {
                        paste0(x,
                            collapse = ","
                        )
                    }
                ))]
            # if user specifies add a column to notify of the imputation
            if (imputation_ind) {
                sumstats_dt[, IMPUTATION_A1 := TRUE]
                sumstats_dt[, IMPUTATION_A2 := TRUE]
            }
        }
        # remove rows where A1/A2 couldn't be found
        # If user wants log, save it to there
        if (log_folder_ind) {
            name <- "alleles_not_found_from_snp"
            name <- get_unique_name_log_file(name = name,
                                             log_files = log_files)
            write_sumstats(
                sumstats_dt =
                    sumstats_dt[!complete.cases(sumstats_dt[, c(
                        "A1",
                        "A2"
                    )]), ],
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
                paste0(check_save_out$log_folder, "/",
                       name, check_save_out$extension)
        }
        sumstats_dt <- sumstats_dt[complete.cases(
            sumstats_dt[, c("A1", "A2")]), ]
        # move SNP, CHR, BP, A1 and A2 to start
        other_cols <-
            names(sumstats_dt)[!names(sumstats_dt) %in%
                c("SNP", "CHR", "BP", "A1", "A2")]
        data.table::setcolorder(
            sumstats_dt,
            c("SNP", "CHR", "BP", "A1", "A2", other_cols)
        )

        return(list(
            "sumstats_dt" = sumstats_dt,
            "rsids" = rsids,
            "allele_flip_check" = allele_flip_check,
            "log_files" = log_files,
            "bi_allelic_filter" = bi_allelic_filter
        ))
    } else {
        return(list(
            "sumstats_dt" = sumstats_dt, 
            "rsids" = rsids,
            "allele_flip_check" = allele_flip_check,
            "log_files" = log_files,
            "bi_allelic_filter" = bi_allelic_filter
        ))
    }
}
