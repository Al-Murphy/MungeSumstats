#' Get VCF sample ID(s)
#'
#' @inheritParams format_sumstats
#' @importFrom VariantAnnotation scanVcfHeader
#' @return sample_id
#' @keywords internal
get_vcf_sample_ids <- function(path) {
    if (startsWith(path, "https://gwas.mrcieu.ac.uk")) {
        header <- VariantAnnotation::scanVcfHeader(file = path)
        sample_id <- header@samples
    } else {
        header <- readLines(path, 100)
        sample_id <- header[grepl("^##SAMPLE", header)] # gets ##SAMPLE
        if (length(sample_id) == 0) {
            ### Try again with more rows
            header <- readLines(path, 500)
            sample_id <- header[grepl("^##SAMPLE", header)] # gets ##SAMPLE
        }
        # get rid of everything after ID
        sample_id <- gsub(",.*$", "", sample_id)
        # rmv ##SAMPLE=
        sample_id <- base::substr(sample_id, 10, nchar(sample_id))
        # rmv things before equals
        sample_id <- sub(".+=(.+)", "\\1", sample_id)
    }
    if (length(sample_id) == 0) {
        message(
            "No ##SAMPLE row found.",
            "Will infer sample name(s) from data colnames."
        )
        return(NULL)
    } else {
        return(sample_id)
    }
}
