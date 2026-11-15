\name{harmonise_file}
\alias{harmonise_file}
\title{Stream a preservation-oriented harmonisation edition to disk}
\usage{
harmonise_file(input, reference, path,
  columns = list(chr = "CHR", pos = "BP", id = "SNP", other = "A1",
    effect_allele = "A2", effect = "BETA", eaf = "FRQ"),
  build = NULL, frequency_type = c("unknown", "EAF", "MAF"),
  effect_type = c("beta", "log_odds", "none"),
  delimiter = "\t", batch_size = 100000L, candidate_limit = 10000L)
}
\arguments{
\item{input}{Plain or gzip input table.}
\item{reference}{An \code{mss_reference} or reference directory.}
\item{path}{New immutable output directory.}
\item{columns,build,frequency_type,effect_type,candidate_limit}{See \code{\link{harmonise_sumstats}}.}
\item{delimiter}{A single tab or comma. Unquoted one-record-per-line tables only.}
\item{batch_size}{Maximum source lines resident in a batch.}
}
\value{A manifest containing output path, source count, per-status audit counts, options and checksums.}
\description{Processes a file in bounded line batches and publishes an output edition only after complete parsing and audit validation. Plain and gzip-compressed inputs are supported.}
\details{The original bytes are retained in \code{source-artifact}; decompressed source lines remain in \code{source.tsv} or \code{source.csv}. Normalized TSV rows carry a 1-based \code{source_row}; numeric fields are serialized with 17 significant digits. Audit counts and an RDS manifest identify every row and the reference version. Portable quoted TSV receipts include \code{metadata.tsv} (key/value), \code{columns.tsv} (field/source_column), \code{reference_checksums.tsv} (name/md5), and \code{checksums.tsv} (filename/md5). The checksum list excludes itself and the RDS manifest. \code{normalized_build} labels the actual native coordinate build. Reference hashes are verified once before streaming. Field-count errors, quotes and read warnings (including missing final newline and detected gzip truncation) abort without publishing a final directory. Syntactically valid incomplete files require an external acquisition receipt to detect. Memory scales with batch size, line lengths and the native candidate cap; neither complete source nor complete output tables are loaded in R. No full-reference production-scale benchmark or Windows execution is claimed by this implementation.}
\seealso{\code{\link{prepare_reference}}, \code{\link{harmonise_sumstats}}}
