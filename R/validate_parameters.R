#' Ensure that the input parameters are logical
#'
#' @inheritParams format_sumstats
#' @return No return
#' @keywords internal
validate_parameters <- function(path,
                                ref_genome,
                                convert_ref_genome,
                                convert_small_p,
                                es_is_beta,
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
                                chr_style,
                                rmv_chr,
                                on_ref_genome,
                                infer_eff_direction,
                                eff_on_minor_alleles,
                                strand_ambig_filter,
                                allele_flip_check,
                                allele_flip_drop,
                                allele_flip_z,
                                allele_flip_frq,
                                bi_allelic_filter,
                                flip_frq_as_biallelic,
                                snp_ids_are_rs_ids,
                                remove_multi_rs_snp,
                                frq_is_maf,
                                indels,
                                drop_indels,
                                check_dups,
                                dbSNP,
                                dbSNP_tarball,
                                write_vcf,
                                return_format,
                                ldsc_format,
                                save_format,
                                imputation_ind,
                                log_folder_ind,
                                log_mungesumstats_msgs,
                                mapping_file,
                                tabix_index,
                                chain_source,
                                local_chain,
                                drop_na_cols,
                                #deprecated parameters
                                rmv_chrPrefix) {
    # Checking if the file exists should happen first -
    # can pass dt/df of sumstats
    pth_msg <- paste0(
        "Path to GWAS sumstats is not valid, pass a file path or a ",
        "dataframe/data.table object to the path parameter"
    )
    if (!is.data.frame(path) && !file.exists(path) &&
        !startsWith(path, "https://gwas.mrcieu.ac.uk") &&
        !startsWith(path, "https://")) {
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
        "The chosen genome build to convert to must be one of ",
        "GRCh37 or GRCh38 ",
        "or left as null so the genome build will be inferred ",
        "from the data."
    )
    if (!is.null(convert_ref_genome) && !(toupper(convert_ref_genome) %in%
        c("GRCH37", "GRCH38"))) {
        stop(gen_msg2)
    }
    #check chain file source
    chain_msg <- paste0(
      "The chosen chain file source to convert to must be one of ",
      "Ensembl or UCSC ('ensembl','ucsc')"
    )
    if(length(chain_source)>1 || !tolower(chain_source) %in% c("ucsc",
                                                               "ensembl")){
      stop(chain_msg)
    }
    #check local chain file
    if (!is.null(local_chain) && !file.exists(local_chain)){
      lcl_chain_msg <- paste0(
        "The local_chain parameter is invalid, please chose a valid path to a ",
        "local chain file or leave as NULL to download a chain file."
      )
      stop(lcl_chain_msg)
    }
    # if no tarball, we require that the matching SNPlocs pkg is installed
    if (is.null(dbSNP_tarball)) {
          if (!is.numeric(dbSNP) || dbSNP <= 0){
            stp_msg <- paste0(
              "`dbSNP` must be a positive integer (e.g. 144, 155...",
              ") or provide `dbSNP_tarball`.")
            stop(stp_msg)
          }
          pkg_name <- sprintf(
                "SNPlocs.Hsapiens.dbSNP%d.GRCh%s",
                as.integer(dbSNP),
                ifelse(toupper(ref_genome)=="GRCH37","37","38")
            )
          if (!requireNamespace(pkg_name, quietly=TRUE)) {
                stop("To use dbSNP=", dbSNP, " on ", ref_genome,
                     ", please install the Bioconductor package '", 
                     pkg_name, "'.", call. = FALSE)
            }
    } else{
      if(!file.exists(dbSNP_tarball)){
        tarball_msg <- paste0(
          "The dbSNP_tarball parameter is invalid, please chose a valid path ",
          "to a local dbSNP release in tarball format or leave as NULL to use ",
          "the dbSNP parameter"
        )
        stop(tarball_msg)
      }
    }
    # checks for installed packages
    GRCH37_msg1 <- paste0(
        "Install 'SNPlocs.Hsapiens.dbSNP144.GRCh37' to use ",
        "'GRCh37' as 'ref_genome' and 144 as dbSNP version"
    )
    GRCH37_msg2 <- paste0(
        "Install 'BSgenome.Hsapiens.1000genomes.hs37d5' to ",
        "use 'GRCh37' as 'ref_genome'"
    )
    GRCH37_msg3 <- paste0(
      "Install 'SNPlocs.Hsapiens.dbSNP155.GRCh37' to ",
      "use 'GRCh37' as 'ref_genome' and 155 as dbSNP version"
    )
    #dbsnp 144 grch37
    if (any(toupper(ref_genome) == "GRCH37") && as.integer(dbSNP)==144 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37",
            quietly = TRUE
        )) {
        stop(GRCH37_msg1)
    }
    #dbsnp 155 grch37
    if (any(toupper(ref_genome) == "GRCH37") && as.integer(dbSNP)==155 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh37",
                          quietly = TRUE
        )) {
      stop(GRCH37_msg3)
    }
    #grch37 ref genome
    if (any(toupper(ref_genome) == "GRCH37") &&
        !requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
            quietly = TRUE
        )) {
        stop(GRCH37_msg2)
    }
    GRCH38_msg1 <- paste0(
        "Install 'SNPlocs.Hsapiens.dbSNP144.GRCh38' to use ",
        "'GRCh38' as 'ref_genome' and 144 as dbSNP version"
    )
    GRCH38_msg2 <- paste0(
        "Install 'BSgenome.Hsapiens.NCBI.GRCh38' to ",
        "use 'GRCh38' as 'ref_genome'"
    )
    GRCH38_msg3 <- paste0(
      "Install 'SNPlocs.Hsapiens.dbSNP155.GRCh38' to ",
      "use 'GRCh38' as 'ref_genome' and 155 as dbSNP version"
    )
    #grch38 and dbsnp 144
    if (any(toupper(ref_genome) == "GRCH38") && as.integer(dbSNP)==144 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38",
            quietly = TRUE
        )) {
        stop(GRCH38_msg1)
    }
    #grch38 and dbsnp 155
    if (any(toupper(ref_genome) == "GRCH38") && as.integer(dbSNP)==155 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh38",
                          quietly = TRUE
        )) {
      stop(GRCH38_msg3)
    }
    #ref genome
    if (any(toupper(ref_genome) == "GRCH38") &&
        !requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38",
            quietly = TRUE
        )) {
        stop(GRCH38_msg2)
    }
    #if ref genome not specified, have to install both
    if (is.null(ref_genome)) {
        #ref genomes
        if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
                            quietly = TRUE
        )) {
          stop(GRCH37_msg2)
        }
        if (!requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38",
                            quietly = TRUE
        )) {
          stop(GRCH38_msg2)
        }
        #dbSNP
        if (as.integer(dbSNP)==144 &&
            !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37",
            quietly = TRUE
        )) {
            stop(GRCH37_msg1)
        }
        if (as.integer(dbSNP)==155 &&
            !requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh37",
                            quietly = TRUE
        )) {
            stop(GRCH37_msg3)
        }
        if (as.integer(dbSNP)==144 &&
            !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38",
            quietly = TRUE
            )) {
            stop(GRCH38_msg1)
        }
        if (as.integer(dbSNP)==155 &&
            !requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh38",
                            quietly = TRUE
            )) {
            stop(GRCH38_msg3)
        }
    }


    # Check binary values
    if (!is.logical(N_dropNA)) {
        stop("N_dropNA must be either TRUE or FALSE")
    }
    if (!is.logical(convert_small_p)) {
        stop("convert_small_p must be either TRUE or FALSE")
    }
    if (!is.logical(es_is_beta)){
      stop("es_is_beta must be either TRUE or FALSE")
    }
    if (!is.logical(compute_z)) {
        if(!toupper(compute_z) %in% c('P','BETA'))
          stop("compute_z must be either TRUE or FALSE")
    }
    if (!is.logical(convert_n_int)) {
        stop("convert_n_int must be either TRUE or FALSE")
    }
    if (!is.logical(on_ref_genome)) {
        stop("on_ref_genome must be either TRUE or FALSE")
    }
    if(!is.logical(infer_eff_direction)){
      stop("infer_eff_direction must be either TRUE or FALSE")
    }
    if(!is.logical(eff_on_minor_alleles)){
      stop("eff_on_minor_alleles must be either TRUE or FALSE")
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
    if (!is.logical(flip_frq_as_biallelic)) {
      stop("flip_frq_as_biallelic must be either TRUE or FALSE")
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
    if (!is.logical(indels)) {
        stop("indels must be either TRUE or FALSE")
    }
    if (!is.logical(drop_indels)) {
      stop("drop_indels must be either TRUE or FALSE")
    }
    if (!is.logical(check_dups)) {
      stop("check_dups must be either TRUE or FALSE")
    }
    if (!is.logical(write_vcf)) {
        stop("write_vcf must be either TRUE or FALSE")
    }
    #LDSC format has been deprecated - use save_format
    if (!isFALSE(ldsc_format)) {
        stop("`ldsc_format` has been deprecated. Use `save_format='LDSC'`.")
    }
    #save_format
    if(!is.null(save_format) &&
       !tolower(save_format) %in% c("ldsc","opengwas")){
      stop("save_format must be NULL or one of LDSC or openGWAS")
    }
    opengws_err <- paste0("IEU OpenGWAS format only available when saving as ",
                          "VCF. Set `write_vcf=True` and rerun ",
                          "`format_sumstats()`")
    if(!is.null(save_format) &&
        tolower(save_format)=="opengwas" & isFALSE(write_vcf))
      stop(opengws_err)
    if (!is.logical(pos_se)) {
        stop("pos_se must be either TRUE or FALSE")
    }
    if (!is.logical(effect_columns_nonzero)) {
        stop("effect_columns_nonzero must be either TRUE or FALSE")
    }
    if (!is.logical(imputation_ind)) {
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
    #SNP level compute_n value - not allowed
    if(length(compute_n)>1 && is.numeric(compute_n)){
      compute_n_msg <- paste0("Multiple integer values for compute_n not ",
                              "supported.\nAdd SNP specific N values before ",
                              "passing to MungeSumstats.")
      stop(compute_n_msg)
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

    # Check that chr_style is a valid choice
    styles <- c("NCBI", "UCSC", "dbSNP", "Ensembl")
    if (!(chr_style %in% styles)) {
      stop("chr_style must be one of ", paste(styles, collapse = ", "))
    }
    # Check that rmv_chr choices are all valid chromosomes
    # according to the Ensembl/NCBI naming style
    chrs <- c(as.character(seq_len(22)), "X", "Y", "MT")
    if (!is.null(rmv_chr)) {
      if (!all(rmv_chr %in% chrs)) {
        stop("rmv_chr choices must be one or more of: \n",
             paste(chrs, collapse = ", "))
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
           !requireNamespace("MatrixGenerics", quietly = TRUE)) ){
        pkgs <- c("Rsamtools","MatrixGenerics")
        missing_pkgs <- pkgs[!pkgs %in% rownames(utils::installed.packages())]
        tbx_msg <- paste0(
            "The following packages must be installed when tabix_index=TRUE:\n",
            paste(" -",missing_pkgs,
                  collapse = "\n"))
        stop(tbx_msg)
    }
    
    # validate drop_na_cols
    if (!is.character(drop_na_cols)) {
      if (!is.null(drop_na_cols)) {
        stop(
          "Parameter `drop_na_cols` should be either a character vector of column names, or `NULL`"
        )
      }
    } 
    
    #deprecated parameters
    if (!is.null(rmv_chrPrefix)) { 
      dep_msg <- paste0(
        "The parameter rmv_chrPrefix is now deprecated, please use chr_style ",
        "instead.\nThe default of rmv_chrPrefix = True will give the same ",
        "result as using chr_style = 'Ensembl'."
        )
      stop(dep_msg)
    }
    rmv_chrPrefix = rmv_chrPrefix
}
