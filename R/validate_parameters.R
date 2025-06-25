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
                                # deprecated parameters
                                rmv_chrPrefix) {
  #### 0) Early tarball  dbSNP sanity checks ####
  if (!is.null(dbSNP_tarball) && !file.exists(dbSNP_tarball)) {
    stop("Tarball not found: ", dbSNP_tarball, call. = FALSE)
  }
  if (!is.numeric(dbSNP) || dbSNP <= 0) {
    stop("`dbSNP` must be a positive integer, even when using a tarball.")
  }
  #### 1) path  basic genome‐build checks ####
  if (!is.data.frame(path) && !file.exists(path) &&
      !startsWith(path, "https://gwas.mrcieu.ac.uk") &&
      !startsWith(path, "https://")) {
    stop("Path to GWAS sumstats is not valid; must be a file or data.frame.")
  }
  
  ## if user supplied a tarball, it must exist on disk
  if (!is.null(dbSNP_tarball) && !file.exists(dbSNP_tarball)) {
    stop("Tarball not found: ", dbSNP_tarball, call. = FALSE)
  }
  ## dbSNP parameter must always be a positive integer
  if (!is.numeric(dbSNP) || dbSNP <= 0) {
    stop("`dbSNP` must be a positive integer, even when using a tarball.")
  }
  
  # — only skip the SNPlocs/BSgenome package checks if a tarball is supplied —
  skip_pkg_checks <- !is.null(dbSNP_tarball)
  
  if (!is.null(ref_genome) &&
      !toupper(ref_genome) %in% c("GRCH37", "GRCH38")) {
    stop("ref_genome must be one of 'GRCh37' or 'GRCh38', or NULL to infer.")
  }
  if (!is.null(convert_ref_genome) &&
      !toupper(convert_ref_genome) %in% c("GRCH37", "GRCH38")) {
    stop("convert_ref_genome must be one of 'GRCh37' or 'GRCh38', or NULL.")
  }
  if (length(chain_source) > 1 ||
      !tolower(chain_source) %in% c("ucsc", "ensembl")) {
    stop("chain_source must be either 'ucsc' or 'ensembl'.")
  }
  if (!is.null(local_chain) && !file.exists(local_chain)) {
    stop("local_chain must be NULL or a valid path to a chain file.")
  }
  
  #### 2) SNPlocs & BSgenome package checks (only when no tarball) ####
  if (!skip_pkg_checks)  {
    # ensure dbSNP is a positive integer
    if (!is.numeric(dbSNP) || dbSNP <= 0) {
      stop("`dbSNP` must be a positive integer (e.g. 144,155,156) or provide dbSNP_tarball.")
    }
    
    # per‐build SNPlocs sanity
    if (toupper(ref_genome) == "GRCH37" && dbSNP == 144 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37", quietly =
                          TRUE)) {
      stop("Install 'SNPlocs.Hsapiens.dbSNP144.GRCh37' for GRCh37dbSNP144.")
    }
    if (toupper(ref_genome) == "GRCH37" && dbSNP == 155 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh37", quietly =
                          TRUE)) {
      stop("Install 'SNPlocs.Hsapiens.dbSNP155.GRCh37' for GRCh37dbSNP155.")
    }
    if (toupper(ref_genome) == "GRCH38" && dbSNP == 144 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38", quietly =
                          TRUE)) {
      stop("Install 'SNPlocs.Hsapiens.dbSNP144.GRCh38' for GRCh38dbSNP144.")
    }
    if (toupper(ref_genome) == "GRCH38" && dbSNP == 155 &&
        !requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh38", quietly =
                          TRUE)) {
      stop("Install 'SNPlocs.Hsapiens.dbSNP155.GRCh38' for GRCh38dbSNP155.")
    }
    
    # BSgenome fallback when ref_genome is NULL
    if (is.null(ref_genome)) {
      if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly =
                            TRUE)) {
        stop("Install 'BSgenome.Hsapiens.1000genomes.hs37d5' for GRCh37 support.")
      }
      if (!requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
        stop("Install 'BSgenome.Hsapiens.NCBI.GRCh38' for GRCh38 support.")
      }
    }
  }
  
  #### 3) all the remaining logical/numeric/mapping‐file checks ####
  if (!is.logical(N_dropNA))
    stop("N_dropNA must be TRUE or FALSE")
  if (!is.logical(convert_small_p))
    stop("convert_small_p must be TRUE or FALSE")
  if (!is.logical(es_is_beta))
    stop("es_is_beta must be TRUE or FALSE")
  if (!is.logical(compute_z))
    stop("compute_z must be TRUE or FALSE or 'P'/'BETA'")
  if (!is.logical(convert_n_int))
    stop("convert_n_int must be TRUE or FALSE")
  if (!is.logical(on_ref_genome))
    stop("on_ref_genome must be TRUE or FALSE")
  if (!is.logical(infer_eff_direction))
    stop("infer_eff_direction must be TRUE or FALSE")
  if (!is.logical(eff_on_minor_alleles))
    stop("eff_on_minor_alleles must be TRUE or FALSE")
  if (!is.logical(strand_ambig_filter))
    stop("strand_ambig_filter must be TRUE or FALSE")
  if (!is.logical(allele_flip_check))
    stop("allele_flip_check must be TRUE or FALSE")
  if (!is.logical(allele_flip_drop))
    stop("allele_flip_drop must be TRUE or FALSE")
  if (!is.logical(allele_flip_z))
    stop("allele_flip_z must be TRUE or FALSE")
  if (!is.logical(allele_flip_frq))
    stop("allele_flip_frq must be TRUE or FALSE")
  if (!is.logical(bi_allelic_filter))
    stop("bi_allelic_filter must be TRUE or FALSE")
  if (!is.logical(flip_frq_as_biallelic))
    stop("flip_frq_as_biallelic must be TRUE or FALSE")
  if (!is.logical(snp_ids_are_rs_ids))
    stop("snp_ids_are_rs_ids must be TRUE or FALSE")
  if (!is.logical(remove_multi_rs_snp))
    stop("remove_multi_rs_snp must be TRUE or FALSE")
  if (!is.logical(frq_is_maf))
    stop("frq_is_maf must be TRUE or FALSE")
  if (!is.logical(indels))
    stop("indels must be TRUE or FALSE")
  if (!is.logical(drop_indels))
    stop("drop_indels must be TRUE or FALSE")
  if (!is.logical(check_dups))
    stop("check_dups must be TRUE or FALSE")
  if (!is.logical(write_vcf))
    stop("write_vcf must be TRUE or FALSE")
  if (!isFALSE(ldsc_format))
    stop("`ldsc_format` is deprecated; use save_format='LDSC'")
  if (!is.null(save_format) &&
      !tolower(save_format) %in% c("ldsc", "opengwas")) {
    stop("save_format must be one of NULL, 'LDSC', or 'openGWAS'")
  }
  if (!is.logical(pos_se))
    stop("pos_se must be TRUE or FALSE")
  if (!is.logical(effect_columns_nonzero))
    stop("effect_columns_nonzero must be TRUE or FALSE")
  if (!is.logical(imputation_ind))
    stop("imputation_ind must be TRUE or FALSE")
  if (!is.logical(log_folder_ind))
    stop("log_folder_ind must be TRUE or FALSE")
  if (!is.logical(log_mungesumstats_msgs))
    stop("log_mungesumstats_msgs must be TRUE or FALSE")
  
  if (!is.numeric(INFO_filter))
    stop("INFO_filter must be numeric")
  if (!is.numeric(FRQ_filter))
    stop("FRQ_filter must be numeric")
  if (!is.numeric(N_std))
    stop("N_std must be numeric")
  
  if (length(compute_n) > 1 && is.numeric(compute_n)) {
    stop("Multiple compute_n integers not supported; supply vector before calling.")
  }
  if (!is.numeric(compute_n) || compute_n < 0) {
    if (is.character(compute_n)) {
      methods <- c("ldsc", "giant", "metal", "sum")
      if (!all(tolower(compute_n) %in% methods)) {
        stop("compute_n must be one of 0 or ",
             paste(methods, collapse = ", "))
      }
    } else {
      stop("compute_n must be 0 or a non‐negative integer.")
    }
  }
  
  styles <- c("NCBI", "UCSC", "dbSNP", "Ensembl")
  if (!chr_style %in% styles) {
    stop("chr_style must be one of ", paste(styles, collapse = ", "))
  }
  chrs <- c(as.character(1:22), "X", "Y", "MT")
  if (!is.null(rmv_chr) && !all(rmv_chr %in% chrs)) {
    stop("rmv_chr must be subset of ", paste(chrs, collapse = ", "))
  }
  
  if (!tolower(return_format) %in%
      c("vr",
        "vranges",
        "gr",
        "granges",
        "genomicranges",
        "data.table")) {
    stop("return_format must be one of 'data.table','GRanges','VRanges'")
  }
  
  ## mapping_file sanity
  essential <- c("SNP", "CHR", "BP", "P", "A1", "A2")
  signed    <- c("Z", "OR", "BETA", "LOG_ODDS", "SIGNED_SUMSTAT")
  if (!is.data.frame(mapping_file) || ncol(mapping_file) != 2 ||
      !all(toupper(colnames(mapping_file)) %in% c("UNCORRECTED", "CORRECTED")) ||
      !all(essential %in% mapping_file$Corrected) ||
      !any(signed   %in% mapping_file$Corrected)) {
    stop("mapping_file must be 2-col DF mapping all essential  ≥1 signed columns.")
  }
  
  if (tabix_index &&
      (
        !requireNamespace("Rsamtools", quietly = TRUE) ||
        !requireNamespace("MatrixGenerics", quietly = TRUE)
      )) {
    stop("tabix_index=TRUE requires Rsamtools and MatrixGenerics installed.")
  }
  
  if (!is.null(drop_na_cols) && !is.character(drop_na_cols)) {
    stop("drop_na_cols must be NULL or a character vector of column names.")
  }
  
  if (!is.null(rmv_chrPrefix)) {
    stop("`rmv_chrPrefix` is deprecated; use chr_style instead.")
  }
}
