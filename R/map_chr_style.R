map_chr_style <- function(CHR, style) {
    chr <- as.character(CHR)
    map_full <- function() GenomeInfoDb::mapSeqlevels(chr, style = style)

    contigs <- unique(chr)
    if (length(contigs) == length(chr) || anyNA(contigs) || any(!nzchar(contigs))) {
        return(map_full())
    }

    # Let the original call handle conditions and unusual mapping results.
    candidates <- tryCatch(
        GenomeInfoDb::mapSeqlevels(contigs, style = style,
                                   best.only = FALSE, drop = FALSE),
        warning = function(warning) NULL,
        error = function(error) NULL
    )
    if (!is.matrix(candidates) || nrow(candidates) == 0L) {
        return(map_full())
    }

    if (length(contigs) == 1L) {
        if (nrow(candidates) != 1L) {
            return(map_full())
        }
        mapped <- rep.int(candidates[[1L]], length(chr))
    } else {
        row_contig <- match(chr, contigs)
        best <- 1L
        if (nrow(candidates) > 1L) {
            # GenomeInfoDb scores maps by matched input rows. Repeated labels
            # contribute once per occurrence.
            counts <- tabulate(row_contig, nbins = length(contigs))
            scores <- drop((!is.na(candidates)) %*% counts)
            best <- which(scores == max(scores))
            if (length(best) != 1L) {
                return(map_full())
            }
        }
        mapped <- unname(candidates[best, ])[row_contig]
    }
    names(mapped) <- chr
    mapped
}
