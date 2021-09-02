#' Check that summary statistics from GWAS are in a homogeneous format
#'
#' @return The address for the modified sumstats file or the actual data
#' dependent on user choice. Also, if log files wanted by the user, the return
#' in both above instances are a list.
#'
#' @examples
#' # Pass path to Educational Attainment Okbay sumstat file to a temp directory
#'
#' eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
#'     package = "MungeSumstats"
#' )
#'
#' ## Call uses reference genome as default with more than 2GB of memory,
#' ## which is more than what 32-bit Windows can handle so remove certain checks
#'
#' is_32bit_windows <-
#'     .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
#' if (!is_32bit_windows) {
#'     reformatted <- format_sumstats(
#'         path = eduAttainOkbayPth,
#'         ref_genome = "GRCh37"
#'     )
#' } else {
#'     reformatted <- format_sumstats(
#'         path = eduAttainOkbayPth,
#'         ref_genome = "GRCh37",
#'         on_ref_genome = FALSE,
#'         strand_ambig_filter = FALSE,
#'         bi_allelic_filter = FALSE,
#'         allele_flip_check = FALSE
#'     )
#' }
#' # returned location has the updated summary statistics file
#' @param path Filepath for the summary statistics file to be formatted. A
#' dataframe or datatable of the summary statistics file can also be passed
#' directly to MungeSumstats using the path parameter.
#' @param ref_genome name of the reference genome used for the GWAS ("GRCh37" or
#' "GRCh38"). Argument is case-insensitive. Default is NULL which infers the
#' reference genome from the data.
#' @param convert_ref_genome name of the reference genome to convert to 
#' ("GRCh37" or "GRCh38"). This will only occur if the current genome build does
#' not match. Default is not to convert the genome build (NULL).
#' @param convert_small_p Binary, should p-values < 5e-324 be converted to 0?
#' Small p-values pass the R limit and can cause errors with LDSC/MAGMA and
#' should be converted. Default is TRUE.
#' @param compute_z Whether to compute Z-score column from P. Default is FALSE.
#' **Note** that imputing the Z-score for every SNP will not correct be
#' perfectly correct and may result in a loss of power. This should only be done
#' as a last resort.
#' @param force_new_z When a "Z" column already exists, it will be used by
#' default. To override and compute a new Z-score column from P set
#' \code{force_new_z=TRUE}.
#' @param compute_n Whether to impute N. Default of 0 won't impute, any other
#' integer will be imputed as the N (sample size) for every SNP in the dataset.
#' **Note** that imputing the sample size for every SNP is not correct and
#' should only be done as a last resort. N can also be inputted with "ldsc", 
#' "sum", "giant" or "metal" by passing one of these for this field or a vector 
#' of multiple. Sum and an integer value creates an N column in the output 
#' whereas giant, metal or ldsc create an Neff or effective sample size. If 
#' multiples are passed, the formula used to derive it will be indicated.
#' @param convert_n_int Binary, if N (the number of samples) is not an integer,
#' should this be rounded? Default is TRUE.
#' @param analysis_trait If multiple traits were studied, name of the trait for
#' analysis from the GWAS. Default is NULL.
#' @param INFO_filter numeric The minimum value permissible of the imputation
#' information score (if present in sumstats file). Default 0.9.
#' @param pos_se Binary Should the standard Error (SE) column be checked to
#' ensure it is greater than 0? Those that are, are removed (if present in
#' sumstats file). Default TRUE.
#' @param effect_columns_nonzero Binary should the effect columns in the data
#' BETA,OR (odds ratio),LOG_ODDS,SIGNED_SUMSTAT be checked to ensure no SNP=0.
#' Those that do are removed(if present in sumstats file). Default FALSE.
#' @param N_std numeric The number of standard deviations above the mean a SNP's
#' N is needed to be removed. Default is 5.
#' @param N_dropNA Drop rows where N is missing.Default is TRUE.
#' @param rmv_chr vector or character The chromosomes on which the SNPs should
#' be removed. Use NULL if no filtering necessary. Default is X, Y and
#' mitochondrial.
#' @param rmv_chrPrefix Remove "chr" or "CHR" from chromosome names. Default is
#' TRUE.
#' @param on_ref_genome Binary Should a check take place that all SNPs are on
#' the reference genome by SNP ID. Default is TRUE.
#' @param strand_ambig_filter Binary Should SNPs with strand-ambiguous alleles
#' be removed. Default is FALSE.
#' @param allele_flip_check Binary Should the allele columns be checked against
#' reference genome to infer if flipping is necessary. Default is TRUE.
#' @param allele_flip_drop Binary Should the SNPs for which neither their A1 or
#' A2 base pair values match a reference genome be dropped. Default is TRUE.
#' @param allele_flip_z Binary should the Z-score be flipped along with effect
#' and FRQ columns like Beta? It is assumed to be calculated off the effect size
#' not the P-value and so will be flipped i.e. default TRUE.
#' @param allele_flip_frq Binary should the frequency (FRQ) column be flipped
#' along with effect and z-score columns like Beta? Default TRUE.
#' @param bi_allelic_filter Binary Should non-biallelic SNPs be removed. Default
#' is TRUE.
#' @param snp_ids_are_rs_ids Binary Should the supplied SNP ID's be assumed to
#' be RS IDs. If not, imputation using the SNP ID for other columns like
#' base-pair position or chromosome will not be possible. If set to FALSE, the
#' SNP RS ID will be imputed from the reference genome if possible. Default is
#' TRUE.
#' @param remove_multi_rs_snp Binary Sometimes summary statistics can have
#' multiple RS IDs on one row (i.e. related to one SNP), for example
#' "rs5772025_rs397784053". This can cause an error so by default, the first
#' RS ID will be kept and the rest removed e.g."rs5772025". If you want to just
#' remove these SNPs entirely, set it to TRUE. Default is FALSE.
#' @param sort_coordinates Whether to sort by coordinates.
#' @param nThread Number of threads to use for parallel processes.
#' @param save_path File path to save formatted data. Defaults to
#' \code{tempfile(fileext=".tsv.gz")}.
#' @param write_vcf Whether to write as VCF (TRUE) or tabular file (FALSE).
#' @param tabix_index Index the formatted summary statistics with
#' \href{http://www.htslib.org/doc/tabix.html}{tabix} for fast querying.
#' @param return_data Return \code{data.table}, \code{GRanges} or \code{VRanges}
#' directly to user. Otherwise, return the path to the save data. Default is
#' FALSE.
#' @param return_format If return_data is TRUE. Object type to be returned
#' ("data.table","vranges","granges").
#' @param ldsc_format Binary Ensure that output format meets all requirements
#' to be fed directly into LDSC without the need for additional munging. Default
#' is FALSE
#' @param log_folder_ind Binary Should log files be stored containing all
#' filtered out SNPs (separate file per filter). The data is outputted in the
#' same format specified for the resulting sumstats file. The only exception to
#' this rule is if output is vcf, then log file saved as .tsv.gz. Default is
#' FALSE.
#' @param log_mungesumstats_msgs Binary Should a log be stored containing all
#' messages and errors printed by MungeSumstats in a run. Default is FALSE
#' @param log_folder Filepath to the directory for the log files and the log of
#' MungeSumstats messages to be stored. Default is a temporary directory.
#' @param imputation_ind Binary Should a column be added for each imputation
#' step to show what SNPs have imputed values for differing fields. This
#' includes a field denoting SNP allele flipping (flipped). On the flipped
#' value, this denoted whether the alelles where switched based on
#' MungeSumstats initial choice of A1, A2 from the input column headers and thus
#' may not align with what the creator intended.**Note** these columns will be
#' in the formatted summary statistics returned. Default is FALSE.
#' @param force_new If a formatted file of the same names as \code{save_path}
#' exists, formatting will be skipped and this file will be imported instead
#' (default). Set \code{force_new=TRUE} to override this.
#' @param mapping_file MungeSumstats has a pre-defined column-name mapping file
#' which should cover the most common column headers and their interpretations.
#' However, if a column header that is in youf file is missing of the mapping we
#' give is incorrect you can supply your own mapping file. Must be a 2 column
#' dataframe with column names "Uncorrected" and "Corrected". See
#' data(sumstatsColHeaders) for default mapping and necessary format.
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setcolorder
#' @importFrom data.table as.data.table
#' @importFrom utils read.table
#' @importFrom utils data
#' @export
format_sumstats <- function(path,
                            ref_genome = NULL,
                            convert_ref_genome = NULL,
                            convert_small_p = TRUE,
                            compute_z = FALSE,
                            force_new_z = FALSE,
                            compute_n = 0L,
                            convert_n_int = TRUE,
                            analysis_trait = NULL,
                            INFO_filter = 0.9,
                            pos_se = TRUE,
                            effect_columns_nonzero = FALSE,
                            N_std = 5,
                            N_dropNA = TRUE,
                            rmv_chr = c("X", "Y", "MT"),
                            rmv_chrPrefix = TRUE,
                            on_ref_genome = TRUE,
                            strand_ambig_filter = FALSE,
                            allele_flip_check = TRUE,
                            allele_flip_drop = TRUE,
                            allele_flip_z = TRUE,
                            allele_flip_frq = TRUE,
                            bi_allelic_filter = TRUE,
                            snp_ids_are_rs_ids = TRUE,
                            remove_multi_rs_snp = FALSE,
                            sort_coordinates = TRUE,
                            nThread = 1,
                            save_path = tempfile(fileext = ".tsv.gz"),
                            write_vcf = FALSE,
                            tabix_index = FALSE,
                            return_data = FALSE,
                            return_format = "data.table",
                            ldsc_format = FALSE,
                            log_folder_ind = FALSE,
                            log_mungesumstats_msgs = FALSE,
                            log_folder = tempdir(),
                            imputation_ind = FALSE,
                            force_new = FALSE,
                            mapping_file = sumstatsColHeaders) {
    #### Setup multi-threading ####
    data.table::setDTthreads(threads = nThread)
    #### Setup empty variables ####
    rsids <- NULL
    orig_dims <- NULL
    log_files <- vector(mode = "list")

    #### Check 1: Ensure save_path is correct.   ####
    check_save_out <- check_save_path(
        save_path = save_path,
        log_folder = log_folder,
        log_folder_ind = log_folder_ind,
        write_vcf = write_vcf
    )
    if (tabix_index && sort_coordinates == FALSE) {
        message(
            "Setting `sort_coordinates=TRUE` in ",
            "order to tabix-index results."
        )
        sort_coordinates <- TRUE
    }

    #### Recognize previously formatted files ####
    if (file.exists(check_save_out$save_path) && force_new == FALSE) {
        message(
            "Importing previously formatted file.",
            "Set `force_new=TRUE` to override this."
        )
        message("    ", check_save_out$save_path)
    } else {
        # Avoid reloading ref genome every time,
        # save it to this parent environment
        # after being made once - speed up code

        # Check input parameters
        validate_parameters(
            path = path,
            ref_genome = ref_genome,
            convert_ref_genome = convert_ref_genome,
            convert_small_p = convert_small_p,
            compute_z = compute_z,
            compute_n = compute_n,
            convert_n_int = convert_n_int,
            analysis_trait = analysis_trait,
            INFO_filter = INFO_filter,
            pos_se = pos_se,
            effect_columns_nonzero = effect_columns_nonzero,
            N_std = N_std,
            N_dropNA = N_dropNA,
            rmv_chr = rmv_chr,
            on_ref_genome = on_ref_genome,
            strand_ambig_filter = strand_ambig_filter,
            allele_flip_check = allele_flip_check,
            allele_flip_drop = allele_flip_drop,
            allele_flip_z = allele_flip_z,
            allele_flip_frq = allele_flip_frq,
            bi_allelic_filter = bi_allelic_filter,
            snp_ids_are_rs_ids = snp_ids_are_rs_ids,
            write_vcf = write_vcf,
            return_format = return_format,
            ldsc_format = ldsc_format,
            imputation_ind = imputation_ind,
            log_folder_ind = log_folder_ind,
            log_mungesumstats_msgs = log_mungesumstats_msgs,
            mapping_file = mapping_file
        )

        # save messages to file if user specified
        if (log_mungesumstats_msgs) {
            msgcon <-
                file(paste0(
                    check_save_out$log_folder,
                    "/MungeSumstats_log_msg.txt"
                ),
                open = "a"
                )
            sink(
                file = paste0(
                    check_save_out$log_folder,
                    "/MungeSumstats_log_output.txt"
                ),
                split = TRUE, append = TRUE
            )
            sink(msgcon, type = "message") # does not support split
            # add name to log_file list
            log_files[["MungeSumstats_log_msg"]] <-
                paste0(
                    check_save_out$log_folder,
                    "/MungeSumstats_log_msg.txt"
                )
            log_files[["MungeSumstats_log_output"]] <-
                paste0(
                    check_save_out$log_folder,
                    "/MungeSumstats_log_output.txt"
                )
        }
        # This almost surely modifies the file (since most sumstats
        # from different studies are differently formatted),
        # so it makes more sense to just make a
        # temporary file <tmp>, and return the address of the temp

        ####  Check 2: Check input format and import ####
        sumstats_return <- list()
        # if data.frame/data.table read it in directly, otherwise read from path
        if (is.data.frame(path)) {
            sumstats_return[["sumstats_dt"]] <- data.table::as.data.table(path)
            # update path in case it causes issue later and for space
            path <- ""
        } else {
            sumstats_return[["sumstats_dt"]] <- read_sumstats(
                path = path,
                nThread = nThread
            )
        }

        #### Check 3:Standardise headers for all OS ####
        sumstats_return <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_return$sumstats_dt,
                mapping_file = mapping_file
            )

        #### If ldsc_format=TRUE, make sure all arguments comply with with.
        check_ldsc <- check_ldsc_format(
            sumstats_dt = sumstats_return$sumstats_dt,
            ldsc_format = ldsc_format,
            convert_n_int = convert_n_int,
            allele_flip_check = allele_flip_check,
            compute_z = compute_z,
            compute_n = compute_n
        )
        convert_n_int <- check_ldsc$convert_n_int
        allele_flip_check <- check_ldsc$allele_flip_check
        compute_z <- check_ldsc$compute_z

        ### Report the number of SNP/CHR/etc. before any filtering
        ### (but after header formatting)
        report_summary(sumstats_dt = sumstats_return$sumstats_dt)
        orig_dims <- dim(sumstats_return$sumstats_dt)

        #### Check 4: Check if multi models used
        # or multi traits tested in GWAS ####
        sumstats_return <-
            check_multi_gwas(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                analysis_trait = analysis_trait,
                mapping_file = mapping_file
            )

        #### Check 33: Check if multi RS ID SNPs in one line ####
        sumstats_return <-
            check_multi_rs_snp(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                remove_multi_rs_snp = remove_multi_rs_snp,
                imputation_ind = imputation_ind,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Infer reference genome if necessary ####
        if (is.null(ref_genome)) {
            ref_genome <- get_genome_build(
                sumstats = sumstats_return$sumstats_dt,
                standardise_headers = FALSE, ## done prev
                sampled_snps = 10000,
                mapping_file = mapping_file
            )
        }

        #### Check 5: Check for uniformity in SNP col - ####
        #### no mix of rs/missing rs/chr:bp ####
        sumstats_return <-
            check_no_rs_snp(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                ref_genome = ref_genome,
                snp_ids_are_rs_ids = snp_ids_are_rs_ids,
                imputation_ind = imputation_ind,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 6: Check for combined allele column (A1 and A2) ####
        sumstats_return <-
            check_allele_merge(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path
            )

        col_headers <- names(sumstats_return$sumstats_dt)

        # Series of checks if CHR or BP columns aren't present
        if (sum(c("CHR", "BP") %in% col_headers) != 2) {
            msg <-
                paste0(
                    "Summary statistics file does not have",
                    "obvious CHR/BP columns.",
                    "Checking to see if they are joined in another column."
                )
            message(msg)

            #### Check 6: check if CHR:BP:A2:A1 merged to 1 column
            sumstats_return <- check_four_step_col(
                sumstats_dt =
                    sumstats_return$sumstats_dt,
                path = path
            )

            #### Check 7: check if there is a column of
            # data with CHR:BP format ####
            sumstats_return <- check_two_step_col(
                sumstats_dt =
                    sumstats_return$sumstats_dt,
                path = path
            )
            #### Re-standardise in case the joined column
            # headers were unusual ####
            sumstats_return <-
                standardise_sumstats_column_headers_crossplatform(
                    sumstats_dt = sumstats_return$sumstats_dt,
                    mapping_file = mapping_file
                )
        }

        #### Check 8: check if CHR and BP are missing but SNP is present ####
        sumstats_return <-
            check_no_chr_bp(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                ref_genome = ref_genome,
                rsids = rsids,
                imputation_ind = imputation_ind,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files
        rsids <- sumstats_return$rsids # update rsids
        sumstats_return$rsids <- NULL

        #### Check 9: check if CHR and BP are present but SNP is missing ####
        sumstats_return <- check_no_snp(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path,
            ref_genome = ref_genome,
            imputation_ind = imputation_ind,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 25: check that all snps are present on reference genome ####
        sumstats_return <- check_on_ref_genome(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            path = path,
            ref_genome = ref_genome,
            on_ref_genome = on_ref_genome,
            rsids = rsids,
            imputation_ind = imputation_ind,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files
        rsids <- sumstats_return$rsids # update rsids
        sumstats_return$rsids <- NULL

        #### Check 10: check if SNP is present but A1 and/or A2 is missing ####
        sumstats_return <-
            check_no_allele(
                sumstats_dt = sumstats_return$sumstats_dt, path = path,
                ref_genome = ref_genome, rsids = rsids,
                imputation_ind = imputation_ind,
                allele_flip_check = allele_flip_check,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files,
                bi_allelic_filter = bi_allelic_filter
            )
        # update values
        log_files <- sumstats_return$log_files
        bi_allelic_filter <- sumstats_return$bi_allelic_filter
        rsids <- sumstats_return$rsids # update rsids
        sumstats_return$rsids <- NULL
        # get updated flip
        allele_flip_check <- sumstats_return$allele_flip_check

        #### Check 11: check that all the vital columns are present ###
        check_vital_col(sumstats_dt = sumstats_return$sumstats_dt)

        #### Check 12: check there is at least one signed sumstats column ###
        check_signed_col(sumstats_dt = sumstats_return$sumstats_dt)

        #### Check 13: check for allele flipping ####
        sumstats_return <-
            check_allele_flip(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                ref_genome = ref_genome,
                rsids = rsids,
                allele_flip_check = allele_flip_check,
                allele_flip_drop = allele_flip_drop,
                allele_flip_z = allele_flip_z,
                allele_flip_frq = allele_flip_frq,
                bi_allelic_filter = bi_allelic_filter,
                imputation_ind = imputation_ind,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files,
                mapping_file = mapping_file
            )
        # update values
        log_files <- sumstats_return$log_files
        rsids <- sumstats_return$rsids # update rsids
        sumstats_return$rsids <- NULL

        #### Check 14: check first three column headers are SNP, CHR, BP ###
        ### (in that order) and also check A1 and A2 are fourth and fifth####
        sumstats_return <-
            check_col_order(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path
            )

        #### Check 15: Keep only rows which have the number
        # of columns expected ####
        sumstats_return <-
            check_miss_data(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 16: check for duplicated columns ####
        # The formatting process can (rarely) result in duplicated columns,
        # i.e. CHR, if CHR:BP is expanded and
        # one already exists...delete duplicates
        sumstats_return <- check_dup_col(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path
        )

        #### Check 17: check for small P-values (3e-400 or lower) ####
        sumstats_return <-
            check_small_p_val(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                convert_small_p = convert_small_p,
                imputation_ind = imputation_ind
            )

        #### Check 18: check is N column not all integers,
        # if so round it up ####
        sumstats_return <-
            check_n_int(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                convert_n_int = convert_n_int,
                imputation_ind = imputation_ind
            )

        #### Check 19: check all rows have SNPs starting with SNP or rs, ####
        #### drop those that don't ####
        sumstats_return <- check_row_snp(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 20: check all rows for duplicated SNPs,
        # remove any that are ####
        sumstats_return <- check_dup_snp(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 21: check all rows for duplicated BPs,
        # remove any that are ####
        sumstats_return <- check_dup_bp(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 22: check for low INFO scores ####
        sumstats_return <-
            check_info_score(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                INFO_filter = INFO_filter,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 30: check standard error is positive ####
        sumstats_return <-
            check_pos_se(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                pos_se = pos_se,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 31: check effect columns are not 0 ####
        sumstats_return <-
            check_effect_columns_nonzero(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                effect_columns_nonzero = effect_columns_nonzero,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 23: check for N > X std dev above mean ####
        sumstats_return <- check_n_num(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path,
            N_std = N_std,
            N_dropNA = N_dropNA,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 24: check that no snps are on specific chromosomes ####
        sumstats_return <- check_chr(
            sumstats_dt = sumstats_return$sumstats_dt,
            path = path,
            rmv_chr = rmv_chr,
            rmv_chrPrefix = rmv_chrPrefix,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 26: check that all snps are not strand ambiguous ####
        sumstats_return <- check_strand_ambiguous(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            path = path,
            ref_genome = ref_genome,
            strand_ambig_filter =
                strand_ambig_filter,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 27: check for non-biallelic SNPS ####
        sumstats_return <- check_bi_allelic(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            path = path,
            ref_genome = ref_genome,
            bi_allelic_filter = bi_allelic_filter,
            rsids = rsids,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files
        )
        # update values
        log_files <- sumstats_return$log_files
        rsids <- sumstats_return$rsids # update rsids
        sumstats_return$rsids <- NULL

        #### Check 28: Compute Z-score ####
        sumstats_return <- check_zscore(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            compute_z = compute_z,
            force_new_z = force_new_z,
            imputation_ind = imputation_ind,
            mapping_file = mapping_file
        )

        #### Check 32: Compute N ####
        sumstats_return <- compute_nsize(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            compute_n = compute_n,
            imputation_ind = imputation_ind
        )

        #### Check 34: Perform liftover ####
        sumstats_return$sumstats_dt <- liftover(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            convert_ref_genome = convert_ref_genome,
            ref_genome = ref_genome,
            imputation_ind = imputation_ind
        )
        
        #### Check 29: Sort rows by genomic coordinates ####
        sumstats_return$sumstats_dt <- sort_coords(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            sort_coordinates =
                sort_coordinates
        )


        #### WRITE data.table TO PATH ####
        write_sumstats(
            sumstats_dt = sumstats_return$sumstats_dt,
            save_path = check_save_out$save_path,
            sep = check_save_out$sep,
            write_vcf = write_vcf,
            tabix_index = tabix_index,
            nThread = nThread
        )
        rm(rsids) # free up memory

        #### Report summary ####
        report_summary(
            sumstats_dt = sumstats_return$sumstats_dt,
            orig_dims = orig_dims
        )
    }



    #### Preview sumstats ####
    preview_sumstats(
        save_path = check_save_out$save_path,
        nrows = 5L
    )

    # if user wanted log of messages remember to unsink at end
    if (log_mungesumstats_msgs) {
        sink(NULL, type = "message")
        sink(NULL, type = "output")
    }

    if (return_data) {
        message("Returning data directly.")
        #### Load data into memory when a pre-existing file is being used
        if (!exists("sumstats_return")) {
            sumstats_return <- list()
            sumstats_return[["sumstats_dt"]] <-
                read_sumstats(
                    path = check_save_out$save_path,
                    nThread = nThread
                )
        }
        out <- convert_sumstats(
            sumstats_dt = sumstats_return$sumstats_dt,
            return_format = return_format
        )
        # if user wants log files return a list
        if (log_folder_ind || log_mungesumstats_msgs) {
            return(list("sumstats" = out, "log_files" = log_files))
        }
        return(out)
    } else {
        message("Returning path to saved data.")
        if (log_folder_ind || log_mungesumstats_msgs) {
            return(list(
                "sumstats" = check_save_out$save_path,
                "log_files" = log_files
            ))
        }
        return(check_save_out$save_path) # Returns address of modified file
    }
}
