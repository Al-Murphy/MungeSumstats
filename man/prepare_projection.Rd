\name{prepare_projection}
\alias{prepare_projection}
\title{Prepare an explicit build projection}
\usage{
prepare_projection(chain, source_reference, target_reference,
  source_build, target_build, source_contigs = character(),
  target_contigs = character())
}
\arguments{
\item{chain}{Plain or gzip UCSC chain with source on its t side and target on its q side.}
\item{source_reference,target_reference}{Prepared portable references.}
\item{source_build,target_build}{Explicit labels that must match the respective references.}
\item{source_contigs,target_contigs}{Optional uniquely named character vectors mapping chain contig names to exact reference names, e.g. \code{c(chr22="22")}. No naming conventions are inferred.}
}
\value{A process-local \code{mss_projection} handle and receipt with checksums, explicit builds, reference versions, mappings and separate algorithm/version identities. Reopen after serialization.}
\description{Loads and validates an explicit chain block index for use by the preservation APIs. The complete native normalized REF span must map to one ungapped block.}
\details{Source and target references are verified; known contig sizes must agree with chain headers. Multiple overlapping mappings, chain gaps, missing mappings, unresolved native identity and target REF conflict retain the source row with missing projected fields. Negative mappings reverse-complement both alleles and validate/left-normalize with the trusted target reference contract. The implementation never assumes the other allele is REF and never infers a target REF/ALT swap. Signed effects retain native ALT direction. Source REF-effect EAF complements additionally require target biallelic catalog evidence; sequence-only or multiallelic target sites cannot supply it. Unchanged source ALT EAF remains unchanged. The native chain block index is resident; variants remain batched. Immutable reference and chain handles verify file size/mtime after initial checksum validation. Mutations preserving both size and mtime are outside this contract. Native fields are retained; separate projected coordinates, alleles, ID, effect, EAF, build, chain ID, strand and status are added. GRCh38 coordinate aliases are populated only for an explicit GRCh38 target. Gap-spanning realignment, reference-switch rescue and specialist transformations remain unsupported and nullable.}
\seealso{\code{\link{harmonise_sumstats}}, \code{\link{harmonise_file}}}
