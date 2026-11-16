#' Prepare an explicit build projection
#'
#' @param chain Plain or gzip UCSC chain, with source on its t side and target on
#'   its q side. No build or contig naming conventions are inferred.
#' @param source_reference,target_reference Prepared portable references.
#' @param source_build,target_build Explicit build labels matching the references.
#' @param source_contigs,target_contigs Optional named character vectors mapping
#'   chain contig names to exact reference contig names, e.g. `c(chr22="22")`.
#' @return A process-local `mss_projection` handle accepted by the preservation APIs.
#' @details Projects the complete native normalized REF span only within one
#'   ungapped chain block. Multiple overlapping chains, gaps, unmapped sequence,
#'   unresolved native identities and target REF conflicts remain nullable.
#'   Negative chains reverse-complement both alleles, then the trusted target
#'   normalization contract validates REF and left aligns. Target REF/ALT swaps
#'   are not inferred. Native ALT effect direction and nullable EAF are retained;
#'   no second EAF complement is attempted. Reference and chain identities are
#'   checksummed; files must stay immutable. The chain block index is resident in
#'   native memory; input variants remain batched. Reopen after serialization.
#' @export
prepare_projection <- function(chain, source_reference, target_reference,
                               source_build, target_build,
                               source_contigs = character(), target_contigs = character()) {
    source <- open_reference(source_reference)
    target <- open_reference(target_reference)
    if (missing(source_build) || missing(target_build) ||
        !identical(source_build, source$build) || !identical(target_build, target$build))
        stop("Explicit source/target build labels must match references")
    validate_names <- function(x) {
        if (!is.character(x) || anyNA(x) || any(!nzchar(x)) ||
            (length(x) && (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
                          anyDuplicated(names(x)) || anyDuplicated(x)))) stop("Contig mappings must be uniquely named character vectors")
    }
    validate_names(source_contigs); validate_names(target_contigs)
    chain <- normalizePath(chain, mustWork = TRUE)
    checksum <- unname(tools::md5sum(chain))
    signature <- file.info(chain)[, c("size", "mtime")]
    reader <- gzfile(chain, "rb"); on.exit(close(reader), add = TRUE)
    plain <- tempfile("mss-chain-"); on.exit(unlink(plain), add = TRUE)
    writer <- file(plain, "wb"); on.exit(close(writer), add = TRUE)
    repeat {
        bytes <- withCallingHandlers(readBin(reader, "raw", 65536L),
                                     warning = function(w) stop("Chain compressed read failure: ", conditionMessage(w)))
        if (!length(bytes)) break
        writeBin(bytes, writer)
    }
    flush(writer)
    contig_sizes <- function(reference) {
        index <- utils::read.delim(file.path(reference$path, "fasta.fai"), header = FALSE,
                                  colClasses = c("character", rep("numeric", 4)))
        stats::setNames(index[[2L]], index[[1L]])
    }
    native <- mss_projection_open_cpp(plain, source_contigs, target_contigs,
                                      contig_sizes(source), contig_sizes(target))
    if (!identical(checksum, unname(tools::md5sum(chain))) ||
        !identical(signature, file.info(chain)[, c("size", "mtime")])) stop("Chain changed during preparation")
    receipt <- list(implementation_version = "mss-projection-2", normalization_implementation = "mss-preservation-7",
                    reference_format = if (identical(source$core_version, target$core_version)) source$core_version else "mixed", chain = chain, chain_md5 = checksum,
                    source_build = source_build, target_build = target_build,
                    source_reference_version = source$version, target_reference_version = target$version,
                    source_checksums = source$checksums, target_checksums = target$checksums,
                    source_contigs = source_contigs, target_contigs = target_contigs,
                    algorithm = "single_ungapped_block_ref_span;target_trusted_ref;no_ref_alt_swap")
    if (source$core_version == "mss-preservation-3" || target$core_version == "mss-preservation-3") {
        receipt$source_reference_format <- source$core_version
        receipt$target_reference_format <- target$core_version
    }
    structure(list(.native = native, source = source, target = target, chain = chain,
                   .signature = signature, receipt = receipt), class = "mss_projection")
}

.mss_project <- function(normalized, reference, projection, candidate_limit) {
    if (!inherits(projection, "mss_projection")) stop("projection must be prepared with prepare_projection")
    if (!identical(reference$checksums, projection$source$checksums) ||
        !identical(reference$build, projection$receipt$source_build)) stop("Projection source reference/build mismatch")
    if (!identical(file.info(projection$chain)[, c("size", "mtime")], projection$.signature)) stop("Projection chain changed; reopen immutable inputs")
    .mss_verify_reference(projection$source, FALSE)
    .mss_verify_reference(projection$target, FALSE)
    result <- mss_project_cpp(projection$.native, projection$target$.native,
                             normalized$chr, normalized$pos, normalized$ref, normalized$alt,
                             normalized$effect, normalized$eaf, candidate_limit, normalized$swapped)
    result$projected_build <- rep(projection$receipt$target_build, nrow(normalized))
    if (identical(projection$receipt$target_build, "GRCh38")) {
        for (field in c("chr", "pos", "ref", "alt")) normalized[[paste0("grch38_", field)]] <- result[[paste0("projected_", field)]]
    }
    cbind(normalized, result)
}
