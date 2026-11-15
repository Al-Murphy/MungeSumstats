\name{harmonise_file}
\alias{harmonise_file}
\title{Stream a preservation-oriented harmonisation edition to disk}
\usage{
harmonise_file(input, reference, path,
  columns = list(chr = "CHR", pos = "BP", id = "SNP", other = "A1",
    effect_allele = "A2", effect = "BETA", eaf = "FRQ"),
  build = NULL, frequency_type = c("unknown", "EAF", "MAF"),
  effect_type = c("beta", "log_odds", "none"),
  delimiter = "\t", batch_size = 100000L, candidate_limit = 10000L,
  source_strand = NULL, projection = NULL, batch_bytes = 8L * 1024L * 1024L,
  max_record_bytes = 64L * 1024L * 1024L, acquisition_receipt = NULL)
}
\arguments{
\item{input}{Plain or gzip input table.}
\item{reference}{An \code{mss_reference} or reference directory.}
\item{path}{New immutable output directory.}
\item{columns,build,frequency_type,effect_type,candidate_limit,source_strand}{See \code{\link{harmonise_sumstats}}. A static strand declaration is recorded in receipts and does not add or modify source columns.}
\item{delimiter}{A single tab or comma. Double-quoted fields, doubled escaped quotes, embedded delimiters/newlines, CRLF and final records without a newline are supported.}
\item{projection}{Optional handle from \code{\link{prepare_projection}}.}
\item{batch_bytes}{Decoded-token batch byte budget; a single larger record is allowed up to max_record_bytes.}
\item{max_record_bytes}{Maximum encoded bytes per record; oversized records abort without publication.}
\item{acquisition_receipt}{Optional named list with expected rows, original file bytes and/or md5. Mismatches abort publication.}
\item{batch_size}{Maximum source records resident in a batch.}
}
\value{A manifest containing output path, source count, per-status audit counts, options and checksums.}
\description{Processes a file in adaptive bounded byte/row batches and publishes an output edition only after complete parsing and audit validation. Plain and gzip-compressed inputs are supported.}
\details{The original bytes are retained in \code{source-artifact}; exact decompressed source bytes remain in \code{source.tsv} or \code{source.csv}. Normalized TSV rows carry a 1-based \code{source_row}; numeric fields are serialized with 17 significant digits. Audit counts and an RDS manifest identify every row and the reference version. Portable quoted TSV receipts include \code{metadata.tsv} (key/value), \code{columns.tsv} (field/source_column), \code{reference_checksums.tsv} (name/md5), and \code{checksums.tsv} (filename/md5). The checksum list excludes itself and the RDS manifest. \code{normalized_build} labels the actual native coordinate build. Reference hashes are verified once before streaming. Malformed quotes, field-count errors and read warnings (including detected gzip truncation) abort without publishing a final directory. Missing final LF is valid syntax. Syntactically valid incomplete files require an external acquisition receipt to detect. The original bytes and parsed token spelling are preserved; quoting syntax is decoded for harmonisation only. Memory scales with batch byte/row limits, individual record limit, parser chunk (at most 64 KiB), per-field R object overhead and native candidate cap; neither complete source nor complete output tables are loaded in R. Projection receipts and a separate projection status audit are included when projection is requested. Projection implementation, normalization implementation and reference format identities remain distinct.}
\seealso{\code{\link{prepare_reference}}, \code{\link{harmonise_sumstats}}}
