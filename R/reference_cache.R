# Active only during format_sumstats(); nested calls restore the previous loader.
.reference_cache <- new.env(parent = emptyenv())
.reference_cache$loader <- NULL

# Cache immutable reference lookup results within one R process.
# Every caller receives a copy because downstream code modifies data.tables.
make_reference_cache <- function(loader) {
  partitions <- new.env(parent=emptyenv())
  function(snps, ref_genome, dbSNP=c(144,155), msg=NULL, chr_filt=NULL) {
    # Unusual IDs and chromosome subsets keep the original lookup semantics.
    supported <- is.character(snps) && !anyNA(snps) &&
      all(grepl('^rs[1-9][0-9]*$', snps)) && is.null(chr_filt) &&
      length(dbSNP)==1L && dbSNP %in% c(144,155) &&
      length(ref_genome)==1L && toupper(ref_genome) %in% c('GRCH37','GRCH38')
    if (!isTRUE(supported)) return(loader(snps=snps, ref_genome=ref_genome, dbSNP=dbSNP, msg=msg, chr_filt=chr_filt))
    key <- paste(toupper(ref_genome), dbSNP, sep=':')
    ids <- unique(snps)
    state <- if (exists(key, envir=partitions, inherits=FALSE)) get(key, envir=partitions) else
      list(queried=character(), data=NULL)
    missing <- setdiff(ids, state$queried)
    fresh <- NULL
    if (length(missing) || is.null(state$data)) {
      fresh <- loader(snps=missing, ref_genome=ref_genome, dbSNP=dbSNP, msg=msg, chr_filt=chr_filt)
      if (!data.table::is.data.table(fresh) || !'SNP' %in% names(fresh) ||
          anyNA(fresh$SNP) || !all(fresh$SNP %in% missing)) {
        # Alias expansion cannot safely use the requested-ID cache key.
        return(loader(snps=snps, ref_genome=ref_genome, dbSNP=dbSNP, msg=msg, chr_filt=chr_filt))
      }
      state$data <- if (is.null(state$data)) data.table::copy(fresh) else
        data.table::rbindlist(list(state$data, fresh), use.names=TRUE)
      data.table::setkeyv(state$data, 'SNP')
      state$queried <- union(state$queried, missing)
      assign(key, state, envir=partitions)
    }
    result <- data.table::copy(state$data[state$data$SNP %in% ids])
    data.table::setkeyv(result, 'SNP')
    result
  }
}
