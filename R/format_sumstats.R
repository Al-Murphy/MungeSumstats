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
#' ## Using dbSNP = 144 for speed as it's smaller but you should use 155 unless
#' ## you know what you are doing and need 144
#'
#' is_32bit_windows <-
#'     .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
#' if (!is_32bit_windows) {
#'     reformatted <- format_sumstats(
#'         path = eduAttainOkbayPth,
#'         ref_genome = "GRCh37",
#'         dbSNP = 144
#'     )
#' } else {
#'     reformatted <- format_sumstats(
#'         path = eduAttainOkbayPth,
#'         ref_genome = "GRCh37",
#'         on_ref_genome = FALSE,
#'         strand_ambig_filter = FALSE,
#'         bi_allelic_filter = FALSE,
#'         allele_flip_check = FALSE,
#'         dbSNP=144
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
#' @param chain_source source of the chain file to use in liftover, if converting
#' genome build ("ucsc" or "ensembl"). Note that the UCSC chain files require a
#' license for commercial use. The Ensembl chain is used by default ("ensembl").
#' @param local_chain Path to local chain file to use instead of downlaoding.
#' Default of NULL i.e. no local file to use. NOTE if passing a local chain file
#' make sure to specify the path to convert from and to the correct build like 
#' GRCh37 to GRCh38. We can not sense check this for local files. The chain file
#' can be submitted as a gz file (as downloaed from source) or unzipped.
#' @param convert_small_p Binary, should non-negative
#' p-values <= 5e-324 be converted to 0?
#' Small p-values pass the R limit and can cause errors with LDSC/MAGMA and
#' should be converted. Default is TRUE.
#' @param convert_large_p Binary, should p-values >1 be converted to 1?
#' P-values >1 should not be possible and can cause errors with LDSC/MAGMA and
#' should be converted. Default is TRUE.
#' @param convert_neg_p Binary, should p-values <0 be converted to 0?
#' Negative p-values should not be possible and can cause errors
#' with LDSC/MAGMA and should be converted. Default is TRUE.
#' @param compute_z Whether to compute Z-score column. Default is FALSE. This
#' can be computed from Beta and SE with (Beta/SE) or P
#' (Z:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))).
#' **Note** that imputing the Z-score from P for every SNP will not be
#' perfectly correct and may result in a loss of power. This should only be done
#' as a last resort. Use 'BETA' to impute by BETA/SE and 'P' to impute by SNP
#' p-value.
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
#' @param impute_beta Binary, whether BETA should be imputed using other effect
#' data if it isn't present in the sumstats. Note that this imputation is an
#' approximation (for Z & SE approach) so could have an effect on downstream
#' analysis. Use with caution. The different methods MungeSumstats will try and
#' impute beta (in this order or priority) are:
#' 1. log(OR)  2. Z x SE
#' Default value is FALSE.
#' @param es_is_beta Binary, whether to map ES to BETA. We take BETA to be any
#' BETA-like value (including Effect Size). If this is not the case for your
#' sumstats, change this to FALSE. Default is TRUE.
#' @param impute_se Binary, whether the standard error should be imputed using
#' other effect data if it isn't present in the sumstats. Note that this
#' imputation is an approximation so could have an effect on downstream
#' analysis. Use with caution. The different methods MungeSumstats will try and
#' impute se (in this order or priority) are:
#' 1. BETA / Z  2. abs(BETA/ qnorm(P/2))
#' Default is FALSE.
#' @param analysis_trait If multiple traits were studied, name of the trait for
#' analysis from the GWAS. Default is NULL.
#' @param ignore_multi_trait If you have multiple traits (p-values) in the study
#' but you want to ignorwe these and instead use a standard named p-value, set
#' to TRUE. By default is FALSE which will check for multi-traits.
#' @param INFO_filter numeric The minimum value permissible of the imputation
#' information score (if present in sumstats file). Default 0.9.
#' @param FRQ_filter numeric The minimum value permissible of the frequency(FRQ)
#' of the SNP (i.e. Allele Frequency (AF)) (if present in sumstats file). By
#' default no filtering is done, i.e. value of 0.
#' @param pos_se Binary Should the standard Error (SE) column be checked to
#' ensure it is greater than 0? Those that are, are removed (if present in
#' sumstats file). Default TRUE.
#' @param effect_columns_nonzero Binary should the effect columns in the data
#' BETA,OR (odds ratio),LOG_ODDS,SIGNED_SUMSTAT be checked to ensure no SNP=0.
#' Those that do are removed(if present in sumstats file). Default FALSE.
#' @param N_std numeric The number of standard deviations above the mean a SNP's
#' N is needed to be removed. Default is 5.
#' @param N_dropNA Drop rows where N is missing.Default is TRUE.
#' @param chr_style Chromosome naming style to use in the formatted summary
#'   statistics file ("NCBI", "UCSC", "dbSNP", or "Ensembl"). The NCBI and
#'   Ensembl styles both code chromosomes as `1-22, X, Y, MT`; the UCSC style is
#'   `chr1-chr22, chrX, chrY, chrM`; and the dbSNP style is
#'   `ch1-ch22, chX, chY, chMT`. Default is Ensembl.
#' @param rmv_chrPrefix Is now deprecated, do. not use. Use chr_style instead -
#' chr_style = 'Ensembl' will give the same result as rmv_chrPrefix=TRUE used to
#' give.
#' @param rmv_chr Chromosomes to exclude from the formatted summary statistics
#'   file. Use NULL if no filtering is necessary. Default is `c("X", "Y", "MT")`
#'   which removes all non-autosomal SNPs.
#' @param on_ref_genome Binary Should a check take place that all SNPs are on
#' the reference genome by SNP ID. Default is TRUE.
#' @param infer_eff_direction Binary Should a check take place to ensure the 
#' alleles match the effect direction? Default is TRUE.
#' @param eff_on_minor_alleles Binary Should MungeSumstats assume that the 
#' effects are majoritively measured on the minor alleles? Default is FALSE as
#' this is an assumption that won't be appropriate in all cases. However, the 
#' benefit is that if we know the majority of SNPs have their effects based on 
#' the minor alleles, we can catch cases where the allele columns have been 
#' mislabelled.
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
#' @param bi_allelic_filter Binary Should non-bi-allelic SNPs be removed.
#' Default is TRUE.
#' @param flip_frq_as_biallelic Binary Should non-bi-allelic SNPs frequency 
#' values be flipped as 1-p despite there being other alternative alleles? 
#' Default is FALSE but if set to TRUE, this allows non-bi-allelic SNPs to be 
#' kept despite needing flipping.
#' @param snp_ids_are_rs_ids Binary Should the supplied SNP ID's be assumed to
#' be RSIDs. If not, imputation using the SNP ID for other columns like
#' base-pair position or chromosome will not be possible. If set to FALSE, the
#' SNP RS ID will be imputed from the reference genome if possible. Default is
#' TRUE.
#' @param remove_multi_rs_snp Binary Sometimes summary statistics can have
#' multiple RSIDs on one row (i.e. related to one SNP), for example
#' "rs5772025_rs397784053". This can cause an error so by default, the first
#' RS ID will be kept and the rest removed e.g."rs5772025". If you want to just
#' remove these SNPs entirely, set it to TRUE. Default is FALSE.
#' @param frq_is_maf Conventionally the FRQ column is intended to show the
#' minor/effect allele frequency (MAF) but sometimes the major allele frequency
#' can be inferred as the FRQ column. This logical variable indicates that the
#' FRQ column should be renamed to MAJOR_ALLELE_FRQ if the frequency values
#' appear to relate to the major allele i.e. >0.5. By default this mapping won't
#' occur i.e. is TRUE.
#' @param indels Binary does your Sumstats file contain Indels? These don't
#' exist in our reference file so they will be excluded from checks if this
#' value is TRUE. Default is TRUE.
#' @param drop_indels Binary, should any indels found in the sumstats be
#' dropped? These can not be checked against a reference dataset and will have
#' the same RS ID and position as SNPs which can affect downstream analysis.
#' Default is False.
#' @param drop_na_cols A character vector of column names to be checked for 
#' missing values. Rows with missing values in any of these columns (if present 
#' in the dataset) will be dropped. If `NULL`, all columns will be checked for 
#' missing values. Default columns are SNP, chromosome, position, allele 1, 
#' allele2, effect columns (frequency, beta, Z-score, standard error, log odds,
#' signed sumstats, odds ratio), p value and N columns. 
#' @param dbSNP version of dbSNP to be used for imputation (144 or 155).
#' @param check_dups whether to check for duplicates - if formatting QTL
#' datasets this should be set to FALSE otherwise keep as TRUE. Default is TRUE.
#' @param sort_coordinates Whether to sort by coordinates of resulting sumstats
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
#' @param ldsc_format DEPRECATED, do not use. Use save_format="LDSC" instead.
#' @param save_format Output format of sumstats. Options are NULL - standardised
#' output format from MungeSumstats, LDSC - output format compatible with LDSC
#' and openGWAS - output compatible with openGWAS VCFs. Default is NULL. 
#' **NOTE** - If LDSC format is used, the naming convention of A1 as the 
#' reference (genome build) allele and A2 as the effect allele will be reversed
#' to match LDSC (A1 will now be the effect allele). See more info on this 
#' [here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68). Note that any 
#' effect columns (e.g. Z) will be inrelation to A1 now instead of A2.
#' @param log_folder_ind Binary Should log files be stored containing all
#' filtered out SNPs (separate file per filter). The data is outputted in the
#' same format specified for the resulting sumstats file. The only exception to
#' this rule is if output is vcf, then log file saved as .tsv.gz. Default is
#' FALSE.
#' @param log_mungesumstats_msgs Binary Should a log be stored containing all
#' messages and errors printed by MungeSumstats in a run. Default is FALSE
#' @param log_folder Filepath to the directory for the log files and the log of
#' MungeSumstats messages to be stored. Default is a temporary directory. Note
#' the name of the log files (log messages and log outputs) are now the same as
#' the name of the file specified in the save path parameter with the extension
#' '_log_msg.txt' and '_log_output.txt' respectively.
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
                            chain_source = "ensembl",
                            local_chain = NULL,
                            convert_small_p = TRUE,
                            convert_large_p = TRUE,
                            convert_neg_p = TRUE,
                            compute_z = FALSE,
                            force_new_z = FALSE,
                            compute_n = 0L,
                            convert_n_int = TRUE,
                            impute_beta = FALSE,
                            es_is_beta = TRUE,
                            impute_se = FALSE,
                            analysis_trait = NULL,
                            ignore_multi_trait = FALSE,
                            INFO_filter = 0.9,
                            FRQ_filter = 0,
                            pos_se = TRUE,
                            effect_columns_nonzero = FALSE,
                            N_std = 5,
                            N_dropNA = TRUE,
                            chr_style = "Ensembl",
                            rmv_chr = c("X", "Y", "MT"),
                            on_ref_genome = TRUE,
                            infer_eff_direction = TRUE,
                            eff_on_minor_alleles = FALSE,
                            strand_ambig_filter = FALSE,
                            allele_flip_check = TRUE,
                            allele_flip_drop = TRUE,
                            allele_flip_z = TRUE,
                            allele_flip_frq = TRUE,
                            bi_allelic_filter = TRUE,
                            flip_frq_as_biallelic = FALSE,
                            snp_ids_are_rs_ids = TRUE,
                            remove_multi_rs_snp = FALSE,
                            frq_is_maf = TRUE,
                            indels = TRUE,
                            drop_indels  = FALSE,
                            drop_na_cols = c("SNP", "CHR", "BP", "A1", "A2", 
                                             "FRQ", "BETA", "Z", "OR", 
                                             "LOG_ODDS", "SIGNED_SUMSTAT", "SE", 
                                             "P", "N"),
                            dbSNP = 155,
                            check_dups = TRUE,
                            sort_coordinates = TRUE,
                            nThread = 1,
                            save_path = tempfile(fileext = ".tsv.gz"),
                            write_vcf = FALSE,
                            tabix_index = FALSE,
                            return_data = FALSE,
                            return_format = "data.table",
                            ldsc_format = FALSE,
                            save_format = NULL,
                            log_folder_ind = FALSE,
                            log_mungesumstats_msgs = FALSE,
                            log_folder = tempdir(),
                            imputation_ind = FALSE,
                            force_new = FALSE,
                            mapping_file = sumstatsColHeaders,
                            #deprecated parameters
                            rmv_chrPrefix = NULL
                            ) {
    #### Setup multi-threading ####
    data.table::setDTthreads(threads = nThread)
    #### Setup empty variables ####
    rsids <- orig_dims <- A1_n <- A2 <- A1 <- NULL
    log_files <- vector(mode = "list")
    t1 <- Sys.time()

    #### Check 1: Ensure save_path is correct.   ####
    check_save_out <- check_save_path(
        save_path = save_path,
        log_folder = log_folder,
        log_folder_ind = log_folder_ind,
        tabix_index = tabix_index,
        write_vcf = write_vcf
    )
    if (isTRUE(tabix_index) && (sort_coordinates == FALSE)) {
        message(
            "Setting `sort_coordinates=TRUE` in ",
            "order to tabix-index results."
        )
        sort_coordinates <- TRUE
    }

    #### Recognize previously formatted files ####
    if (file.exists(check_save_out$save_path) && (force_new == FALSE)) {
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
            es_is_beta = es_is_beta,
            compute_z = compute_z,
            compute_n = compute_n,
            convert_n_int = convert_n_int,
            analysis_trait = analysis_trait,
            INFO_filter = INFO_filter,
            FRQ_filter = FRQ_filter,
            pos_se = pos_se,
            effect_columns_nonzero = effect_columns_nonzero,
            N_std = N_std,
            N_dropNA = N_dropNA,
            chr_style = chr_style,
            rmv_chr = rmv_chr,
            on_ref_genome = on_ref_genome,
            infer_eff_direction = infer_eff_direction,
            eff_on_minor_alleles = eff_on_minor_alleles,
            strand_ambig_filter = strand_ambig_filter,
            allele_flip_check = allele_flip_check,
            allele_flip_drop = allele_flip_drop,
            allele_flip_z = allele_flip_z,
            allele_flip_frq = allele_flip_frq,
            bi_allelic_filter = bi_allelic_filter,
            flip_frq_as_biallelic = flip_frq_as_biallelic,
            snp_ids_are_rs_ids = snp_ids_are_rs_ids,
            remove_multi_rs_snp = remove_multi_rs_snp,
            frq_is_maf = frq_is_maf,
            indels = indels,
            drop_indels = drop_indels,
            dbSNP = dbSNP,
            check_dups = check_dups,
            write_vcf = write_vcf,
            return_format = return_format,
            ldsc_format = ldsc_format,
            save_format = save_format,
            imputation_ind = imputation_ind,
            log_folder_ind = log_folder_ind,
            log_mungesumstats_msgs = log_mungesumstats_msgs,
            mapping_file = mapping_file,
            tabix_index = tabix_index,
            chain_source = chain_source,
            local_chain = local_chain,
            drop_na_cols = drop_na_cols,
            #deprecated parameters
            rmv_chrPrefix = rmv_chrPrefix
        )

        # save messages to file if user specified
        if (log_mungesumstats_msgs) {
            #get name of file from save_path
            nme <- strsplit(basename(check_save_out$save_path),
                            split="[.]")[[1]][1]
            msg <- paste0("Saving output messages to:\n",
                          paste0(check_save_out$log_folder,"/",
                                 nme,"_log_msg.txt"),"\n",
                          "Any runtime errors will be saved to:\n",
                          paste0(check_save_out$log_folder,"/",
                                 nme,"_log_output.txt"),"\n",
                          "Messages will not be printed to terminal.")
            message(msg)
            msgcon <-
                file(paste0(
                  check_save_out$log_folder,"/",
                  nme,"_log_msg.txt"),
                open = "a"
                )
            sink(
                file = paste0(
                  check_save_out$log_folder,"/",
                  nme,"_log_output.txt"
                ),
                split = TRUE, append = TRUE
            )
            sink(msgcon, type = "message") # does not support split
            # add name to log_file list
            log_files[["MungeSumstats_log_msg"]] <-
                paste0(
                  check_save_out$log_folder,"/",
                  nme,"_log_msg.txt"
                )
            log_files[["MungeSumstats_log_output"]] <-
                paste0(
                  check_save_out$log_folder,"/",
                  nme,"_log_output.txt"
                )
        }
        # This almost surely modifies the file (since most sumstats
        # from different studies are differently formatted),
        # so it makes more sense to just make a
        # temporary file <tmp>, and return the address of the temp

        #Ensure dbSNP is a integer (make using it later easier)
        #already validated in validate param function
        dbSNP <- as.integer(dbSNP)

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
                samples = if(is.null(analysis_trait)) 1 else analysis_trait,
                nThread = nThread
            )
        }
        
        #If user inputted mapping file, validate
        if(!identical(mapping_file, sumstatsColHeaders)) {
            message("Non-standard mapping file detected.",
                    "Making sure all entries in `Uncorrected`",
                    " are in upper case.")
            data.table::setDF(mapping_file)
            #check again
            if(!identical(mapping_file, sumstatsColHeaders)) {
              mapping_file$Uncorrected <- toupper(mapping_file$Uncorrected)
            }  
        }
        
        #If es_is_beta remove from mapping file if present
        if (!es_is_beta & nrow(mapping_file[mapping_file$Uncorrected=="ES" &
                                           mapping_file$Corrected=="BETA",])>=1)
          {
          mapping_file <- mapping_file[!(mapping_file$Uncorrected=="ES" &
                                         mapping_file$Corrected=="BETA"),]
          #Add ES mapping
          es_cols <- data.frame("Uncorrected"=c("ES","EFFECT_SIZE",
                                                "EFFECT.SIZE","EFFECT-SIZE",
                                                "EFFECT SIZE",
                                                "EFFECT_SIZE_ESTIMATE",
                                                "EFFECT SIZE ESTIMATE",
                                                "EFFECT.SIZE.ESTIMATE",
                                                "ES.A1","ES.A2","ES-A1","ES-A2",
                                                "ES_A1","ES_A2"),
                                "Corrected"=rep("ES",14))
          mapping_file <- rbind(mapping_file,es_cols)
        }
        
        #### Check 40:Check for log10 p instead of p ####
        sumstats_return <-
          read_log_pval(sumstats_dt = sumstats_return$sumstats_dt)
 
        #### Check 2:Check for effect direction ####
        sumstats_return <-
          infer_effect_column(
            sumstats_dt = sumstats_return$sumstats_dt,
            mapping_file = mapping_file,
            dbSNP = dbSNP,
            nThread = nThread,
            ref_genome = ref_genome,
            on_ref_genome = on_ref_genome,
            infer_eff_direction = infer_eff_direction,
            eff_on_minor_alleles = eff_on_minor_alleles
          )
        
        #### Check 3:Standardise headers for all OS ####
        sumstats_return <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_return$sumstats_dt,
                mapping_file = mapping_file
            )

        #### If save_format=LDSC, make sure all arguments comply with with.
        check_ldsc <- check_ldsc_format(
            sumstats_dt = sumstats_return$sumstats_dt,
            save_format = save_format,
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
                ignore_multi_trait = ignore_multi_trait,
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

        #before running inference of genome build, do any formatting
        #not using the reference sets

        #### Infer reference genome if necessary ####
        if (is.null(ref_genome)) {
            ref_genome <- get_genome_build(
                sumstats = sumstats_return$sumstats_dt,
                standardise_headers = FALSE, ## done prev
                sampled_snps = 10000,
                mapping_file = mapping_file,
                dbSNP=dbSNP
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
                indels=indels,
                imputation_ind = imputation_ind,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files,
                dbSNP = dbSNP
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

        #Ensure A1 and A2 are upper case
        sumstats_return <-
            make_allele_upper(sumstats_dt = sumstats_return$sumstats_dt,
                              log_files = log_files)
        # update values
        log_files <- sumstats_return$log_files

        # Series of checks if CHR or BP columns aren't present
        if (sum(c("CHR", "BP") %in% col_headers) != 2) {
            msg <-
                paste0(
                    "Summary statistics file does not have",
                    " obvious CHR/BP columns. ",
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

        #### Check 38: validate BP
        sumstats_return <- check_bp_range(
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
                log_files = log_files,
                dbSNP = dbSNP
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
            indels = indels,
            imputation_ind = imputation_ind,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files,
            dbSNP = dbSNP
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
            indels=indels,
            rsids = rsids,
            imputation_ind = imputation_ind,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files,
            dbSNP = dbSNP
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
                bi_allelic_filter = bi_allelic_filter,
                dbSNP = dbSNP
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
        sumstats_return <-
            check_signed_col(
                sumstats_dt = sumstats_return$sumstats_dt,
                impute_beta = impute_beta,
                log_folder_ind = log_folder_ind,
                rsids = rsids,
                imputation_ind = imputation_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files
        rsids <- sumstats_return$rsids # update rsids
        sumstats_return$rsids <- NULL

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
                flip_frq_as_biallelic = flip_frq_as_biallelic,
                imputation_ind = imputation_ind,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files,
                mapping_file = mapping_file,
                dbSNP = dbSNP
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
                log_files = log_files,
                drop_na_cols = drop_na_cols
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

        #### Check 17: check for small P-values (<=5e-324) ####
        sumstats_return <-
            check_small_p_val(
                sumstats_dt = sumstats_return$sumstats_dt,
                convert_small_p = convert_small_p,
                imputation_ind = imputation_ind
            )

        #### Check 17.5: check for large (>1) and neg (<0)  p-values ####
        sumstats_return <-
            check_range_p_val(
                sumstats_dt = sumstats_return$sumstats_dt,
                convert_large_p = convert_large_p,
                convert_neg_p = convert_neg_p,
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
        #### drop those that don't ####.
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

        ### Check 37: Drop Indels ###
        sumstats_return <- check_drop_indels(
          sumstats_dt = sumstats_return$sumstats_dt,
          drop_indels=drop_indels,
          path = path,
          log_folder_ind = log_folder_ind,
          check_save_out = check_save_out,
          tabix_index = tabix_index,
          nThread = nThread,
          log_files = log_files)

        # update values
        log_files <- sumstats_return$log_files

        #### Check 20: check all rows for duplicated SNPs,
        # remove any that are ####
        sumstats_return <- check_dup_snp(
            sumstats_dt = sumstats_return$sumstats_dt,
            indels=indels,
            path = path,
            log_folder_ind = log_folder_ind,
            check_save_out = check_save_out,
            tabix_index = tabix_index,
            nThread = nThread,
            log_files = log_files,
            bi_allelic_filter = bi_allelic_filter,
            check_dups = check_dups
        )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 21: check all rows for duplicated BPs,
        # remove any that are ####
        sumstats_return <- check_dup_bp(
            sumstats_dt = sumstats_return$sumstats_dt,
            bi_allelic_filter=bi_allelic_filter,
            check_dups = check_dups,
            indels=indels,
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
                INFO_filter = INFO_filter,
                log_folder_ind = log_folder_ind,
                check_save_out = check_save_out,
                tabix_index = tabix_index,
                nThread = nThread,
                log_files = log_files
            )
        # update values
        log_files <- sumstats_return$log_files

        #### Check 35: check for low FRQ scores ####
        sumstats_return <-
            check_frq(
                sumstats_dt = sumstats_return$sumstats_dt,
                path = path,
                FRQ_filter = FRQ_filter,
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
                impute_se = impute_se,
                imputation_ind = imputation_ind,
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

        #### Check 24: standardize the CHR column ####
        sumstats_return <- check_chr(
          sumstats_dt = sumstats_return$sumstats_dt,
          log_files = log_files,
          check_save_out = check_save_out,
          rmv_chr = rmv_chr,
          nThread = nThread,
          tabix_index = tabix_index,
          log_folder_ind = log_folder_ind
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
            log_files = log_files,
            dbSNP = dbSNP
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

        #### Check 36: Ensure FRQ is MAF ####
        sumstats_return$sumstats_dt <- check_frq_maf(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            frq_is_maf = frq_is_maf
        )


        #### Check 34: Perform liftover ####
        sumstats_return$sumstats_dt <- liftover(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            convert_ref_genome = convert_ref_genome,
            ref_genome = ref_genome,
            imputation_ind = imputation_ind,
            chain_source = chain_source,
            local_chain = local_chain
        )
        #update ref genome of data
        if(!is.null(convert_ref_genome))
          ref_genome <- convert_ref_genome

        #### Check 29: Sort rows by genomic coordinates ####
        sumstats_return$sumstats_dt <- sort_coords(
            sumstats_dt =
                sumstats_return$sumstats_dt,
            sort_coordinates =
                sort_coordinates
        )

        ### Check 39: Ensure CHR follows the requested style ###
        CHR <- NULL
        sumstats_return$sumstats_dt[, 
                                    CHR := GenomeInfoDb::mapSeqlevels(as.character(CHR), 
                                                                      style = chr_style)]
        
        ### IF LDSC, rename A1 and A2, effect columns are fine
        if (!is.null(save_format) && 
            tolower(save_format)=="ldsc") {
          message("Renaming A1,A2 to match LDSC format.")
          #For LDSC format, rename A1 and A2 as LDSC expects A1 to be the effect 
          #column rather than A2 (the opposite to MSS's default) - see more 
          #[here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68).Although, 
          #this didn't seem to make any difference to results in tests, see more
          #https://github.com/neurogenomics/MungeSumstats/issues/160#issuecomment-1891899253
          sumstats_return$sumstats_dt[,A1_n:=A2]
          sumstats_return$sumstats_dt[,A2:=A1]
          sumstats_return$sumstats_dt[,A1:=A1_n]
          sumstats_return$sumstats_dt[,A1_n:=NULL]
        }

        #### WRITE data.table TO PATH ####
        check_save_out$save_path <- write_sumstats(
            sumstats_dt = sumstats_return$sumstats_dt,
            save_path = check_save_out$save_path,
            ref_genome = ref_genome,
            sep = check_save_out$sep,
            write_vcf = write_vcf,
            save_format = save_format,
            tabix_index = tabix_index,
            nThread = nThread,
            return_path = TRUE
        )
        rm(rsids) # free up memory

        #### Report summary ####
        report_summary(
            sumstats_dt = sumstats_return$sumstats_dt,
            orig_dims = orig_dims
        )
    }
    #### Record time taken ####
    t2 <- Sys.time()
    message(
        "Done munging in ",
        round(difftime(t2, t1, units = "mins"), 3), " minutes."
    )

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
            return(list("sumstats" = out,
                        "log_files" = log_files))
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
