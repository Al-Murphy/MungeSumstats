#' Read in VCF format
#'
#' @inheritParams format_sumstats
#' @param temp_save Save temprorary file before proceeding,
#' @param keep_extra_cols Keep non-standard columns.
#' @param nrows integer. The (maximal) number of lines to read.
#' If \code{Inf}, will read in all rows.
#'
#' @return The modified sumstats_file
#' @keywords internal
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
#' @importFrom dplyr %>% rename
#' @importFrom utils head
#' @importFrom stringr str_split
read_vcf <- function(path,
                     nThread = 1,
                     temp_save = FALSE,
                     keep_extra_cols = FALSE,
                     save_path = tempfile(),
                     nrows = Inf) {
    ########## OTHER VCF READERS/WRITERS ###########
    # 1. vcfR
    # 2. VariantAnnotation
    # 3. seqminer
    # 4. Rsamtools
    # 5. echotabix (GH repo: RajLabMSSM/echotabix)
    ################################## 

    ### Add this to avoid confusing BiocCheck
    INFO <- Pval <- P <- LP <- AF <- NULL

    message("Reading VCF file.")
    # fileformat <- gsub("^##fileformat=","",header[1])
    # #First get the name of data column, held in the ##SAMPLE row
    sample_id <- get_vcf_sample_ids(path = path)

    #### Read in full data ####

    #### If the VCF is remotely stored, data.table automatically
    # downloads it to a tmpdir without telling you where.
    # This lets you know where.
    if (startsWith(path, "https://gwas.mrcieu.ac.uk")) {
        vcf_suffixes <- supported_suffixes(
            tabular = FALSE,
            tabular_compressed = FALSE
        )
        tmpdir <-
            file.path(tempdir(), gsub(
                paste(vcf_suffixes, collapse = "|"), "",
                basename(path)
            ), "")
        dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        message("Downloading VCF ==> ", tmpdir)
    } else {
        tmpdir <- tempdir()
    }
    
    #### Read in the data contents of VCF ####
    sumstats_file <- read_vcf_data(path = path, 
                                   nThread = nThread, 
                                   tmpdir = tmpdir,
                                   nrows = nrows) 
    #### Infer sample ID from data colnames if necessary ####
    sample_id <- infer_vcf_sample_ids(
        sample_id = sample_id,
        sumstats_file = sumstats_file
    )
    
    is_parsed <- is_vcf_parsed(sumstats_file = sumstats_file, 
                               verbose = TRUE)
    sumstats_file <- remove_nonstandard_vcf_cols(
        sample_id = sample_id, 
        sumstats_file = sumstats_file,
        keep_extra_cols = keep_extra_cols, 
        standardise_headers = TRUE)
    #### standardise_headers makes sample_id uppercase ####
    sample_id <- toupper(sample_id)
    ## Get format of FORMAT col for parsing
    format <- sumstats_file$FORMAT[1]
    if(!is_parsed && !is.null(format)){   
        format_cols <- stringr::str_split(format, ":")[[1]] 
    }else {
        format_cols <- NULL
    }
    

    # if sample_id col present, split it out into separate columns
    message("Parsing '", sample_id, "' data column.")
    if (length(sample_id) != 0) {
        # split out format into separate values
        # format <- strsplit(format,":")[[1]]
        # First, there is an issue where INFO=".", AF from FORMAT is missing
        # Need to impute 0 for these values
        # See dataset for example of this:
        # http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz
        if ("AF" %in% format &&
            "INFO" %in% names(sumstats_file) &&
            nrow(sumstats_file[INFO == "."]) > 0) {
            # get positions of ":" separators for each SNP
            find_splits <- gregexpr(":", sumstats_file[, get(sample_id)])
            # get position of AF to be imputed
            AF_pos <- which("AF" == format)
            # Update sumstats_file where number of splits isn't correct
            # get row identifiers
            update_rows <-
                (lengths(regmatches(
                    sumstats_file[, get(sample_id)],
                    find_splits
                )) != length(format) - 1)
            # get char pos for imputation
            find_splits <- unlist(lapply(find_splits,
                                         function(x) x[AF_pos - 1]))
            # add to dt
            sumstats_file[, find_splits := find_splits]
            # Now update these by imputing 0 for AF
            sumstats_file[update_rows, (sample_id) :=
                paste0(
                    substr(get(sample_id), 0, find_splits),
                    "0:", # AF impute
                    substr(
                        get(sample_id), find_splits + 1,
                        nchar(get(sample_id))
                    )
                )]
            sumstats_file[, find_splits := NULL]
            # Then replace INFO="." with 0 at later stage
        }
        # check if any cols already present - ID likely will be, rename if so
        if (any(format %in% names(sumstats_file))) {
            format_copy <- format[format %in% names(sumstats_file)]
            for (format_i in format_copy) {
                # If it is ID just remove other ID column
                if (format_i == "ID") {
                    sumstats_file[, (format_i) := NULL]
                } else { # otherwise keep both columns just rename one
                    data.table::setnames(
                        sumstats_file,
                        format_i, paste0(format_i, "2")
                    )
                }
            }
        }
        if(sample_id %in% colnames(sumstats_file) && 
           !is.null(format_cols)){
            sumstats_file[, (format_cols) :=
                              data.table::tstrsplit(get(sample_id),
                                                    split = ":", fixed = TRUE
                              )]
            # Now remove sample_id column
            sumstats_file[, (sample_id) := NULL]
        }
    }
    # sumstatsColHeaders contains mappings for
    # ID to SNP
    # EZ to Z
    # NC to N_CAS
    # SS to N
    # AF to FRQ
    # ES to BETA* -need confirmation

    #### Check for empty columns ####
    empty_cols <- check_empty_cols(
        sumstats_file = sumstats_file,
        sampled_rows = 1000
    )


    if ("MarkerName" %in% colnames(sumstats_file)) {
        if ("ID" %in% empty_cols) {
            message("Replacing empty ID col with MarkerName col.")
            sumstats_file$ID <- sumstats_file$MarkerName
            sumstats_file[, ("MarkerName") := NULL]
        }
    }

    # Need to convert P-value, currently -log10
    if ("Pval" %in% colnames(sumstats_file)) {
        sumstats_file <- sumstats_file %>% dplyr::rename(P = Pval)
    }
    if ("LP" %in% names(sumstats_file)) {
        if ("LP" %in% empty_cols) {
            message("LP column is empty. Cannot compute raw p-values.")
        } else {
            msg <- paste0(
                "VCF file has -log10 P-values, these will be ",
                "converted to unadjusted p-values in the 'P' column."
            )
            message(msg)
            sumstats_file[, P := 10^(-1 * as.numeric(LP))]
        }
    }

    # Need to remove "AF=" at start of INFO column and replace any "." with 0
    if ("INFO" %in% names(sumstats_file)) {
        message("Formatting INFO column.")
        # if info col = "AF=..." then it isn't info, rename to AF if available
        # check first 10k only, or all if less
        num_check <- min(nrow(sumstats_file), 10000)
        # if more than half are, take the column as AF
        if (sum(grepl("^AF=", sumstats_file$INFO)) > num_check / 2) {
            message("INFO column is actually AF, it will be converted.")
            # don't overwirte AF column if it exists
            if (!"AF" %in% names(sumstats_file)) {
                sumstats_file[, AF := INFO]
                sumstats_file[, AF := gsub("^AF=", "", AF)]
            }
        }
        sumstats_file[, INFO := gsub("^AF=", "", INFO)]
        sumstats_file[INFO == ".", INFO := 0]
        # update to numeric
        sumstats_file[, INFO := as.numeric(INFO)]
    }
    if ("INFO" %in% names(empty_cols) && "INFO" %in% names(sumstats_file)) {
        message("NOTE: All INFO scores are empty. Replacing all with 1.")
        sumstats_file$INFO <- 1
    }
    # VCF format has dups of each row, get unique
    sumstats_file <- unique(sumstats_file)
    # #write new data
    if (temp_save) {
        message("Storing intermediate file before proceeding ==> ", path)
        data.table::fwrite(
            x = sumstats_file,
            file = save_path,
            nThread = nThread,
            sep = "\t"
        )
    }
    return(sumstats_file)
}
