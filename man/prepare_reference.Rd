\name{prepare_reference}
\alias{prepare_reference}
\title{Prepare an immutable portable reference index}
\usage{
prepare_reference(fasta, variants, path, build, version, batch_size = 100000L,
                  contig_aliases = character())
}
\arguments{
\item{fasta}{Plain FASTA with fixed-width sequence lines except final lines. Compressed FASTA is unsupported.}
\item{variants}{TSV containing exactly \code{chr,pos,ref,alt,id,aliases}, in that order. Positions are 1-based. Comma-separated ALT alleles and semicolon-separated ID aliases are supported; missing values use a dot.}
\item{path}{New version-specific reference directory; existing paths are refused.}
\item{build}{Explicit coordinate build label.}
\item{version}{Explicit reference release label.}
\item{batch_size}{Maximum index sort-buffer rows; merge fan-in is 32.}
\item{contig_aliases}{Optional uniquely named character vector from authoritative input labels to exact FASTA contig names. Many aliases may name one contig. Unknown destinations, missing-value labels and redirection of existing FASTA names are refused. No labels or assemblies are inferred.}
}
\value{An \code{mss_reference} with paths, build, version and MD5 checksums.}
\description{Builds disk-backed one-to-many coordinate, identifier and alias indexes with a portable C++17 core. Each sequence-resolved REF/ALT pair is trimmed and left aligned against the declared FASTA. The legacy \code{format_sumstats()} interface is unchanged.}
\details{Reference TSV contig names must agree exactly with FASTA. Optional input aliases are stored in checksummed \code{contig_aliases.tsv}; preservation file editions carry \code{reference_contigs.tsv} and the \code{mss-contigs-1} mapping receipt. Only the native-call chromosome labels are mapped; original source fields and bytes remain unchanged. Empty mapping retains strict-name behavior. Source reference tokens remain in \code{source.tsv}; rejected reference alleles are reported in \code{rejected.tsv}. FASTA remains at its declared absolute path and must not be mutated. Index construction uses a bounded row buffer, a merge fan-in of 32 and logarithmically many temporary run paths. Individual source-line lengths are not capped. Full-genome reference-scale memory and performance targets have not been validated.}
\seealso{\code{\link{harmonise_sumstats}}, \code{\link{harmonise_file}}}
