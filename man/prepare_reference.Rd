\name{prepare_reference}
\alias{prepare_reference}
\title{Prepare an immutable portable reference index}
\usage{
prepare_reference(fasta, variants, path, build, version, batch_size = 100000L)
}
\arguments{
\item{fasta}{Plain FASTA with fixed-width sequence lines except final lines. Compressed FASTA is unsupported.}
\item{variants}{TSV containing exactly \code{chr,pos,ref,alt,id,aliases}, in that order. Positions are 1-based. Comma-separated ALT alleles and semicolon-separated ID aliases are supported; missing values use a dot.}
\item{path}{New version-specific reference directory; existing paths are refused.}
\item{build}{Explicit coordinate build label.}
\item{version}{Explicit reference release label.}
\item{batch_size}{Maximum index sort-buffer rows; merge fan-in is 32.}
}
\value{An \code{mss_reference} with paths, build, version and MD5 checksums.}
\description{Builds disk-backed one-to-many coordinate, identifier and alias indexes with a portable C++17 core. Each sequence-resolved REF/ALT pair is trimmed and left aligned against the declared FASTA. The legacy \code{format_sumstats()} interface is unchanged.}
\details{Contig names must agree exactly. Source reference tokens remain in \code{source.tsv}; rejected reference alleles are reported in \code{rejected.tsv}. FASTA remains at its declared absolute path and must not be mutated. Index construction uses a bounded row buffer, a merge fan-in of 32 and logarithmically many temporary run paths. Individual source-line lengths are not capped. Reference-scale memory and performance targets have not been validated.}
\seealso{\code{\link{harmonise_sumstats}}, \code{\link{harmonise_file}}}
