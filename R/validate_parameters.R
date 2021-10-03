#' Ensure that the input parameters are logical
#'
#' @inheritParams format_sumstats
#' @return No return
#' @keywords internal
validate_parameters <- function(path,
                                ref_genome,
                                convert_ref_genome,
                                convert_small_p,
                                compute_z,
                                compute_n,
                                convert_n_int,
                                analysis_trait,
                                INFO_filter,
                                FRQ_filter,
                                pos_se,
                                effect_columns_nonzero,
                                N_std,
                                N_dropNA,
                                rmv_chr,
                                on_ref_genome,
                                strand_ambig_filter,
                                allele_flip_check,
                                allele_flip_drop,
                                allele_flip_z,
                                allele_flip_frq,
                                bi_allelic_filter,
                                snp_ids_are_rs_ids,
                                remove_multi_rs_snp,
                                frq_is_maf,
                                write_vcf,
                                return_format,
                                ldsc_format,
                                imputation_ind,
                                log_folder_ind,
                                log_mungesumstats_msgs,
                                mapping_file,
                                tabix_index) {
    # Checking if the file exists should happen first - 
    # can pass dt/df of sumstats
    pth_msg <- paste0(
        "Path to GWAS sumstats is not valid, pass a file path or a",
        "dataframe/data.table object to the path parameter"
    )
    if (!is.data.frame(path) && !file.exists(path) &&
        !startsWith(path, "https://gwas.mrcieu.ac.uk")) {
        stop(pth_msg)
    }

    gen_msg <- paste0(
        "The chosen genome build must be one of GRCh37 or GRCh38 ",
        "or left as null so the genome build will be inferred ",
        "from the data."
    )
    # Check genome build is valid option
    if (!is.null(ref_genome) && !(toupper(ref_genome) %in%
        c("GRCH37", "GRCH38"))) {
        stop(gen_msg)
    }

    gen_msg2 <- paste0(
        "The chosen genome build to convert to must be one of GRCh37 or GRCh38 ",
        "or left as null so the genome build will be inferred ",
        "from the data."
    )

    if (!is.null(convert_ref_genome) && !(toupper(convert_ref_genome) %in%
        c("GRCH37", "GRCH38"))) {
        stop(gen_msg2)
    }

    # checks for installed packages
    GRCH37_msg1 <- paste0(
        "Install 'SNPlocs.Hsapiens.dbSNP144.GRCh37' to use ",
        "'GRCh37' as 'ref_genome'"
    )
    GRCH37_msg2 <- paste0(
        "Install 'BSgenome.Hsapiens.1000genomes.hs37d5' to ",
        "use 'GRCh37' as 'ref_genome'"
    )
    if (any(toupper(ref_genome) == "GRCH37") &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37",
            quietly = TRUE
        )) {
        stop(GRCH37_msg1)
    }
    if (any(toupper(ref_genome) == "GRCH37") &&
        !requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
            quietly = TRUE
        )) {
        stop(GRCH37_msg2)
    }

    GRCH38_msg1 <- paste0(
        "Install 'SNPlocs.Hsapiens.dbSNP144.GRCh38' to use ",
        "'GRCh38' as 'ref_genome'"
    )
    GRCH38_msg2 <- paste0(
        "Install 'BSgenome.Hsapiens.NCBI.GRCh38' to ",
        "use 'GRCh38' as 'ref_genome'"
    )
    if (any(toupper(ref_genome) == "GRCH38") &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38",
            quietly = TRUE
        )) {
        stop(GRCH38_msg1)
    }
    if (any(toupper(ref_genome) == "GRCH38") &&
        !requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38",
            quietly = TRUE
        )) {
        stop(GRCH38_msg2)
    }

    if (is.null(ref_genome)) {
        if (!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37",
            quietly = TRUE
        )) {
            stop(GRCH37_msg1)
        }
        if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
            quietly = TRUE
        )) {
            stop(GRCH37_msg2)
        }
        if (!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38",
            quietly = TRUE
        )) {
            stop(GRCH38_msg1)
        }
        if (!requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38",
            quietly = TRUE
        )) {
            stop(GRCH38_msg2)
        }
    }


    # Check binary values
    if (!is.logical(N_dropNA)) {
        stop("N_dropNA must be either TRUE or FALSE")
    }
    if (!is.logical(convert_small_p)) {
        stop("convert_small_p must be either TRUE or FALSE")
    }
    if (!is.logical(compute_z)) {
        stop("compute_z must be either TRUE or FALSE")
    }
    if (!is.logical(convert_n_int)) {
        stop("convert_n_int must be either TRUE or FALSE")
    }
    if (!is.logical(on_ref_genome)) {
        stop("on_ref_genome must be either TRUE or FALSE")
    }
    if (!is.logical(strand_ambig_filter)) {
        stop("strand_ambig_filter must be either TRUE or FALSE")
    }
    if (!is.logical(allele_flip_check)) {
        stop("allele_flip_check must be either TRUE or FALSE")
    }
    if (!is.logical(allele_flip_drop)) {
        stop("allele_flip_drop must be either TRUE or FALSE")
    }
    if (!is.logical(allele_flip_z)) {
        stop("allele_flip_z must be either TRUE or FALSE")
    }
    if (!is.logical(allele_flip_frq)) {
        stop("allele_flip_frq must be either TRUE or FALSE")
    }
    if (!is.logical(bi_allelic_filter)) {
        stop("bi_allelic_filter must be either TRUE or FALSE")
    }
    if (!is.logical(snp_ids_are_rs_ids)) {
        stop("snp_ids_are_rs_ids must be either TRUE or FALSE")
    }
    if (!is.logical(remove_multi_rs_snp)) {
        stop("remove_multi_rs_snp must be either TRUE or FALSE")
    }
    if (!is.logical(frq_is_maf)) {
        stop("frq_is_maf must be either TRUE or FALSE")
    }
    if (!is.logical(write_vcf)) {
        stop("write_vcf must be either TRUE or FALSE")
    }
    if (!is.logical(ldsc_format)) {
        stop("ldsc_format must be either TRUE or FALSE")
    }
    if (!is.logical(pos_se)) {
        stop("pos_se must be either TRUE or FALSE")
    }
    if (!is.logical(effect_columns_nonzero)) {
        stop("effect_columns_nonzero must be either TRUE or FALSE")
    }
    if (!is.logical(effect_columns_nonzero)) {
        stop("imputation_ind must be either TRUE or FALSE")
    }
    if (!is.logical(log_folder_ind)) {
        stop("log_folder_ind must be either TRUE or FALSE")
    }
    if (!is.logical(log_mungesumstats_msgs)) {
        stop("log_mungesumstats_msgs must be either TRUE or FALSE")
    }

    # Check numeric
    if (!is.numeric(INFO_filter)) {
        stop("INFO_filter must be numeric")
    }
    if (!is.numeric(FRQ_filter)) {
        stop("FRQ_filter must be numeric")
    }
    if (!is.numeric(N_std)) {
        stop("N_std must be numeric")
    }
    if (!is.numeric(compute_n) || compute_n < 0) {
        if (is.character(compute_n)) {
            methods <- c("ldsc", "giant", "metal", "sum")
            msg_methods <- paste0(
                "Valid methods to compute N are: ",
                paste(methods, collapse = ", "),
                ". Please update compute_n."
            )
            if (!all(tolower(compute_n) %in% methods)) {
                  stop(msg_methods)
              }
        } else {
            stop("compute_n must be 0 or an integer value")
        }
    }

    # Check rmv_chr choices all valid chromosomes
    chrs <- c(as.character(seq_len(22)), "X", "Y", "MT")
    chr_msg <-
        paste0(
            "rmv_chr choices must be one/or more of: \n",
            paste(chrs, collapse = ", ")
        )
    if (!is.null(rmv_chr)) {
        if (!all(rmv_chr %in% chrs)) {
            stop(chr_msg)
        }
    }
    # check return_format
    rf_msg <- paste0(
        "MungeSumstats can only output as data.table, GRanges or ",
        "VRanges objects, not ", return_format,
        ". Check return_format"
    )
    if (!(tolower(return_format) %in% c(
        "vr", "vranges", "gr", "granges",
        "genomicranges", "data.table"
    ))) {
        stop(rf_msg)
    }

    # check if mapping file correct - in case user defined
    essential_cols <- c("SNP", "CHR", "BP", "P", "A1", "A2")
    signed_cols <- c("Z", "OR", "BETA", "LOG_ODDS", "SIGNED_SUMSTAT")
    mapping_file_msg <- paste0(
        "User supplied mapping file should be a 2 column ",
        "dataframe with columns Corrected and Uncorrected",
        ".\nMappings for at least  all essential and one ",
        "signed column must be included.\nThe essential ",
        "columns are: ",
        paste(essential_cols, collapse = ","), ".\n",
        "The signed columns are: ",
        paste(signed_cols, collapse = ","), ".\n",
        "If your inputted column headers ",
        "already have these headers correctly, make sure ",
        "to include a mapping for these too.\nFor example",
        ", SNP -> SNP. See data(sumstatsColHeaders) for ",
        "format."
    )
    if (!is.data.frame(mapping_file) || ncol(mapping_file) != 2 ||
        !all(
            toupper(
                colnames(mapping_file)
            ) %in% c("UNCORRECTED", "CORRECTED")
        ) ||
        !all(essential_cols %in%
            mapping_file[, toupper(colnames(mapping_file)) == "CORRECTED"]) ||
        !any(signed_cols %in%
            mapping_file[, toupper(colnames(mapping_file)) == "CORRECTED"])) {
        stop(mapping_file_msg)
    }
    
    if(tabix_index && 
       any(!requireNamespace("Rsamtools", quietly = TRUE),
           !requireNamespace("seqminer", quietly = TRUE)) ){
        tbx_msg <- paste(
            "Rsamtools and seqminer must be installed when tabix_index=TRUE.")
        stop(tbx_msg)
    }
}
