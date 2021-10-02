#' Ensure A1 & A2 are correctly named, if GWAS SNP constructed as
#' Alternative/Reference or Risk/Nonrisk alleles these SNPs will need to be
#' converted to Reference/Alternative or Nonrisk/Risk. Here non-risk is defined
#' as what's on the reference genome (this may not always be the case).
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @param standardise_headers Run
#' \code{standardise_sumstats_column_headers_crossplatform} first.
#'
#' @return A list containing two data tables:
#' \itemize{
#'   \item \code{sumstats_dt}: the modified summary statistics
#'    \code{data.table} object.
#'   \item \code{rsids}: snpsById, filtered to SNPs of interest if
#'    loaded already. Or else NULL.
#'   \item \code{log_files}: log file list
#' }
#' @keywords internal
#' @import R.utils
#' @importFrom data.table setkey
#' @importFrom data.table :=
#' @importFrom data.table setnames
#' @importFrom data.table set
#' @importFrom data.table setorder
#' @importFrom data.table copy
check_allele_flip <- function(sumstats_dt, path,
                              ref_genome, rsids,
                              allele_flip_check,
                              allele_flip_drop,
                              allele_flip_z,
                              allele_flip_frq,
                              bi_allelic_filter,
                              imputation_ind,
                              log_folder_ind,
                              check_save_out,
                              tabix_index,
                              nThread,
                              log_files,
                              standardise_headers = FALSE,
                              mapping_file) {
    # # GenomicSEM' allele flipping strategy:
    # file.path("https://github.com/GenomicSEM/GenomicSEM/blob",
    #           "fc8f17a817a8022d6900acf41824d27b3676f9c4/R/munge.R#L151")

    # #example
    # path <- system.file("extdata","eduAttainOkbay.txt",
    # package="MungeSumstats")
    # sumstats_dt <- read_sumstats(path = path)
    # sumstats_return <- check_allele_flip(sumstats_dt = sumstats_dt,
    #                                      path=path,
    #                                      ref_genome="GRCh37",
    #                                      rsids=NULL,
    #                                      allele_flip_check=TRUE,
    #                                      standardise_headers=TRUE)
    ## Set variables to be used in in place data.table functions to NULL
    ## to avoid confusing BiocCheck.
    SNP <- i.seqnames <- CHR <- BP <- i.pos <- LP <- P <- A1 <- A2 <- eff_i <-
        i.A1 <- i.A2 <- ss_A1 <- ss_A2 <- i.ref_allele <-
        ref_gen_allele <- match_type <-
        tmp <- flipped <- NULL
    if (standardise_headers) {
        sumstats_dt <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_dt,
                mapping_file = mapping_file
            )[["sumstats_dt"]]
    }
    # If SNP present but no A1/A2 then need to find them
    col_headers <- names(sumstats_dt)
    if (sum(c("A1", "A2") %in% col_headers) == 2 && allele_flip_check) {
        message(
            "Checking for correct direction of A1 (reference) ",
            "and A2 (alternative allele)."
        )
        # check if rsids loaded if not do so
        if (is.null(rsids)) {
            rsids <-
                load_ref_genome_data(
                    data.table::copy(sumstats_dt$SNP),
                    ref_genome, NULL
                )
        }
        # ensure rsids is up-to-date with filtered sumstats_dt
        rsids <- rsids[unique(sumstats_dt$SNP), , nomatch = NULL]
        data.table::setkey(rsids, SNP)
        data.table::setkey(sumstats_dt, SNP)
        # Append reference genome to data and check matches to both A1 or A2
        # For each SNP, if ref genome matches A1, leave it (TRUE)
        # if ref genome matches A2, flip it (FALSE)
        # if ref genome doesn't match A1 or A2, leave it (TRUE)
        sumstats_dt[rsids, ref_gen_allele := i.ref_allele]
        sumstats_dt[is.na(ref_gen_allele), match_type := TRUE]
        sumstats_dt[A1 == ref_gen_allele, match_type := TRUE]
        sumstats_dt[A2 == ref_gen_allele, match_type := FALSE]
        # drop cases that don't match either
        if (allele_flip_drop &&
            nrow(sumstats_dt[A1 != ref_gen_allele &
                A2 != ref_gen_allele, ]) > 0) {
            print_msg0 <-
                paste0(
                    "There are ",
                    formatC(nrow(sumstats_dt[A1 != ref_gen_allele &
                        A2 != ref_gen_allele, ]), big.mark = ","),
                    " SNPs where neither A1 nor A2 match the reference genome.",
                    " These will be removed."
                )
            message(print_msg0)
            # If user wants log, save it to there
            if (log_folder_ind) {
                name <- "alleles_dont_match_ref_gen"
                name <- get_unique_name_log_file(
                    name = name,
                    log_files = log_files
                )
                write_sumstats(
                    sumstats_dt = sumstats_dt[(A1 != ref_gen_allele &
                        A2 != ref_gen_allele), ],
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
            sumstats_dt <- sumstats_dt[!(A1 != ref_gen_allele &
                A2 != ref_gen_allele), ]
        } else {
            sumstats_dt[
                A1 != ref_gen_allele & A2 != ref_gen_allele,
                match_type := TRUE
            ]
        }
        # continue if flipping necessary
        if (nrow(sumstats_dt[match_type == FALSE, ]) > 0) {
            print_msg <- paste0(
                "There are ", 
                formatC(nrow(sumstats_dt[match_type == FALSE, ]),
                        big.mark = ","),
                " SNPs where A1 doesn't match the reference genome.",
                "\nThese will be flipped with their effect columns."
            )
            message(print_msg)
            # swap A1 and A2 for those SNPs needing to be flipped
            sumstats_dt[match_type == FALSE, tmp := A2]
            sumstats_dt[match_type == FALSE, A2 := A1]
            sumstats_dt[match_type == FALSE, A1 := tmp]
            sumstats_dt[, tmp := NULL]

            # flip effect column(s) - BETA, OR, z, log_odds, SIGNED_SUMSTAT, FRQ
            effect_columns <- c("BETA", "OR", "LOG_ODDS", "SIGNED_SUMSTAT")
            if (allele_flip_z) {
                effect_columns <- c(effect_columns, "Z")
            }
            if (allele_flip_frq) {
                effect_columns <- c(effect_columns, "FRQ")
            }
            effect_columns <-
                effect_columns[effect_columns %in% names(sumstats_dt)]
            # if FRQ present and needs to be flipped for SNPs ensure
            # bi_allelic_filter is TRUE
            stp_msg <- paste0(
                "Certain SNPs need to be flipped along with their ",
                "effect columns and frequency column. However to flip ",
                "the\nFRQ column, only bi-allelic SNPs can be ",
                "considered. It is recommended to set ",
                "bi_allelic_filter to TRUE so\nnon-bi-allelic SNPs are",
                " removed. Otherwise, set allele_flip_frq to FALSE to ",
                "not flip the FRQ column but note\nthis could lead to ",
                "incorrect FRQ values."
            )
            if (nrow(sumstats_dt[match_type == FALSE, ]) > 0 &&
                "FRQ" %in% effect_columns &&
                !bi_allelic_filter) {
                stop(stp_msg)
            }
            for (eff_i in effect_columns) { # set updates quicker for DT
                # conversion done in case, VCF beta column may not be numeric
                if (eff_i == "FRQ") {
                    sumstats_dt[
                        match_type == FALSE,
                        (eff_i) := (1 - as.numeric(get(eff_i)))
                    ]
                } else {
                    sumstats_dt[
                        match_type == FALSE,
                        (eff_i) := as.numeric(get(eff_i)) * -1
                    ]
                }
            }
            if (imputation_ind) {
                sumstats_dt[match_type == FALSE, flipped := TRUE]
            }
        }
        # remove extra created columns and return
        sumstats_dt[, ref_gen_allele := NULL]
        sumstats_dt[, match_type := NULL]

        return(list(
            "sumstats_dt" = sumstats_dt,
            "rsids" = rsids, "log_files" = log_files
        ))
    }
    return(list(
        "sumstats_dt" = sumstats_dt,
        "rsids" = rsids, "log_files" = log_files
    ))
}
