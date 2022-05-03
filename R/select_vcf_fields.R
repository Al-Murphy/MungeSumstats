#' Select VCF fields
#' 
#' Select non-empty columns from each VCF field type.
#' 
#' @param nrows Number of rows to use when inferring empty columns.
#' @param verbose Print messages.
#' @inheritParams read_vcf
#' @inheritParams VariantAnnotation::ScanVcfParam
#' 
#' @returns \code{ScanVcfParam} object.
#' @keywords internal
select_vcf_fields <- function(path,
                              nrows=1e7,
                              which=NULL,
                              single_sample=TRUE,
                              nThread=1,
                              verbose=TRUE){
    #### Read header ####
    ## Read the first nrows to determine which columns are useful.
    messager("Finding empty VCF columns based on first",nrows,"rows.",
             v=verbose) 
    #### Make sure file is compressed and indexed ####
    ## File must be indexed in order to use param 
    ## (even if only specifying columns) 
    path <- index_vcf(path = path, 
                      verbose = verbose) 
    has_chr <- infer_chrom_format(path = path, 
                                  nrows = nrows, 
                                  verbose = verbose)
    # #### Construct ScanVcfParam object ####
    param <- VariantAnnotation::ScanVcfParam(
        which = GenomicRanges::GRanges(
            paste0(if(isTRUE(has_chr)) paste0("chr",1) else 1,
                   ":1-",as.integer(nrows))
        )
    )
    vcf_top <- VariantAnnotation::readVcf(file = path,
                                          param = param)
    # #### Warnings generated bc 0 rows are present ####
    df <- suppressWarnings(
        vcf2df(vcf = vcf_top,
               add_sample_names = FALSE)
    )
    #### Check n rows returned ####
    if(nrow(df)==0) {
        stop("Query returned no rows. Increase nrows.")
    }
    #### Check for empty cols #####
    empty_cols <- check_empty_cols(sumstats_dt = df) 
    fields <- VariantAnnotation::vcfFields(x = path)
    for(x in names(fields)){
        fields[[x]] <- fields[[x]][
            !toupper(fields[[x]]) %in% toupper(names(empty_cols))
            ]
    } 
    #### Select samples ####
    samples <- colnames(vcf_top)
    if(length(samples)>1 && isTRUE(single_sample)){
        messager(length(samples),"detected.",
                 "Only using first:",samples[1])
        fields$samples <- fields$samples[fields$samples %in% samples[1]]
    }
    messager("Constructing ScanVcfParam object.")
    #### Select non-empty columns ####
    param <- VariantAnnotation::ScanVcfParam(
        fixed = fields$fixed, 
        info = fields$info, 
        geno = fields$geno, 
        samples = fields$samples,
        trimEmpty = TRUE)
    #### Add genomic ranges if supplied ####
    ## By default, is IRanges::IRangesList() of length 0  
    if(!is.null(which)){
        VariantAnnotation::vcfWhich(param) <- which
    }
    return(param)
}
