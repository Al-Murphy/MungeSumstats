#' Select VCF fields
#' 
#' Select non-empty columns from each VCF field type.
#' @param verbose Print messages.
#' @inheritParams read_vcf
#' @inheritParams check_empty_cols
#' @inheritParams VariantAnnotation::ScanVcfParam
#' 
#' @returns \code{ScanVcfParam} object.
#' @keywords internal
#' @importFrom VariantAnnotation scanVcfHeader ScanVcfParam 
#' @importFrom VariantAnnotation readVcf vcfFields vcfWhich
#' @importFrom GenomicRanges GRanges
#' @importFrom methods is
select_vcf_fields <- function(path,
                              sampled_rows=1e3,
                              which=NULL,
                              samples=NULL,
                              nThread=1,
                              verbose=TRUE){
    #### Read header ####
    ## Read the first n rows to determine which columns are useful.
    messager("Finding empty VCF columns based on first",
             formatC(sampled_rows,big.mark = ","),"rows.",
             v=verbose)  
    #### Get header to extract first chromosome name ####
    # header <- VariantAnnotation::scanVcfHeader(file = path)
    # # #### Construct ScanVcfParam object ####
    # param <- VariantAnnotation::ScanVcfParam(
    #     which = GenomicRanges::GRanges(
    #         paste0(header@reference[[1]],":1-",as.integer(sampled_rows))
    #     )
    # )
    vcf_file <- VariantAnnotation::VcfFile(file = path,
                                           index = paste0(path,".tbi"),
                                           yieldSize=sampled_rows) 
    vcf_top <- VariantAnnotation::readVcf(file = vcf_file)
    fields <- VariantAnnotation::vcfFields(x = vcf_top)
    # #### Warnings generated bc 0 rows are present ####
    df <-  vcf2df(vcf = vcf_top,
                  add_sample_names = FALSE,
                  verbose = verbose)
    remove(vcf_top)
    #### Check n rows returned ####
    if(nrow(df)==0) {
        stop("Query returned no rows. Increase sampled_rows.")
    }
    #### Check for empty cols #####
    remove_empty_cols(sumstats_dt = df)
    for(x in names(fields)[names(fields)!="samples"]){
        fields[[x]] <- fields[[x]][
            fields[[x]] %in% names(df)
            ]
    } 
    #### Select samples ####
    all_samples <- fields$samples
    if(!is.null(samples)){
        #### use first sample only ####
        if(samples==1){
            if(length(all_samples)>1){
                messager(length(all_samples),"samples detected.",
                         "Only using first:",all_samples[1]) 
                all_samples <- all_samples[1] 
            } else if (length(all_samples)==1){
                messager("1 sample detected:",all_samples) 
            } 
        #### User-specified samples #####
        } else if(methods::is(samples,"character")){
            all_samples <- all_samples[
                toupper(all_samples) %in% toupper(samples)
            ]
        }
    }
    ## Filter fields$samples 
    fields$samples <- fields$samples[
        fields$samples %in% all_samples
    ] 
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
