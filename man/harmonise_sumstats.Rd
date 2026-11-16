\name{harmonise_sumstats}
\alias{harmonise_sumstats}
\title{Preserve source rows while harmonising sequence-resolved contrasts}
\usage{
harmonise_sumstats(sumstats, reference,
  columns = list(chr = "CHR", pos = "BP", id = "SNP", other = "A1",
    effect_allele = "A2", effect = "BETA", eaf = "FRQ"),
  build = NULL, frequency_type = c("unknown", "EAF", "MAF"),
  effect_type = c("beta", "log_odds", "none"),
  batch_size = 100000L, candidate_limit = 10000L, verify = TRUE,
  source_strand = NULL, projection = NULL)
}
\arguments{
\item{sumstats}{A data.frame, data.table or previous \code{mss_harmonised} result.}
\item{reference}{An \code{mss_reference} or its directory path.}
\item{columns}{Named list mapping input columns to \code{chr,pos,id,other,effect_allele,effect,eaf}. Missing mappings or columns remain unknown. A2 is the effect allele by default. Optional \code{source_ref,source_alt} mappings declare trusted source REF/ALT, together with either a \code{source_strand} column mapping or the static argument.}
\item{source_strand}{Optional explicit table-level \code{"+"}, \code{"-"} or \code{"unknown"} declaration for trusted source REF/ALT. Cannot be combined with a strand column mapping. REF/ALT and measured alleles must share that strand; positions always identify the forward reference's leftmost 1-based coordinate. A1/A2 do not imply a source-reference declaration.}
\item{projection}{Optional process-local handle from \code{\link{prepare_projection}}. Native fields remain unchanged; separate projected fields and statuses are added.}
\item{build}{Input coordinate build; defaults to reference build. Different builds are rejected.}
\item{frequency_type}{Only explicitly declared EAF is transformed. MAF and unknown frequencies remain missing in normalized output.}
\item{effect_type}{Signed additive beta or log odds, or none. Odds ratios are not accepted.}
\item{batch_size}{Maximum native input batch rows.}
\item{candidate_limit}{Maximum raw candidates per lookup and unique IDs, events and contexts in ambiguity-proof closure. Additional closure work is bounded by 16 times this limit; overflow fails rather than choosing a candidate.}
\item{verify}{Verify all reference MD5 hashes. FALSE only checks sizes and requires already trusted immutable files. Full verification sequentially reads the reference files.}
}
\value{An \code{mss_harmonised} list containing exact \code{source}, row-aligned nullable \code{normalized}, and \code{provenance}.}
\description{A separate preservation interface that retains every input row and original measurement token. It uses a portable native core and disk-indexed reference without calling the legacy filtering pipeline.}
\details{Match, swap, complement and swapped-complement are explicit. Effects are negated only for a uniquely compatible swapped contrast. Frequencies are complemented only when no third allele is known at the reference site; zero is a valid EAF. Palindromic or multi-location ambiguity remains unresolved. ALT-ALT diagnostics require both alleles in one FASTA-valid original REF context, consistent with known input chromosome and position after coherent joint trimming or equivalent indel reanchoring. Shared alternate backgrounds remain part of the contrast even when their strings differ from REF. A competing ALT-ALT explanation also keeps a unique generic REF-ALT hit ambiguous; an identifier does not declare source strand or reference-allele role. ALT-ALT contrasts, symbolic alleles, missing identifiers and malformed optional values remain in source with nullable normalized fields. No arbitrary first ALT is selected. Generic complex alleles require a reported source position even when an identifier names one canonical event; missing-position complex rows return \code{source_coordinates_missing}. Single-base SNP ID inference remains supported. Ambiguity evidence closes all current IDs and exact-event coordinate bridges from the initial coordinate and input-ID pools, retaining each paired original REF context. A distinct proof-only REF-ALT or paired ALT-ALT explanation can withhold a selected mapping but cannot supply a replacement. Unlinked overlapping contexts are outside this explicitly eligible proof universe. Compatible necessary ALT-pair evidence or unknown ALT sequence without a complete background proof defers a selected generic mapping as \code{source_allele_context_required}; this does not assert a proven ALT-ALT contrast. \code{candidate_count} counts the main canonical mapping pool, not proof-graph size or catalog-wide absence. See the installed \code{PRESERVATION_CONTRACT.md} for the full recipe and resource bounds. Projection fields remain missing unless an explicit projection handle is supplied. Optional column mappings \code{entity_type,entity_id} preserve typed specialist entities. Allowed types are \code{small_variant,symbolic_sv,breakend,hla,repeat,gene_burden}. Only symbolic SV and breakend allele syntax is automatically typed. Specialist rows have \code{specialist_preserved} status and nullable small-variant fields; algorithms for specialist transformations are not supplied. Source and output tables are resident in R; use \code{harmonise_file} for a file-to-file bounded-batch path.}
\seealso{\code{\link{prepare_reference}}, \code{\link{harmonise_file}}}
\note{Explicit source REF must validate against FASTA. Source-reference mismatch, unknown strand and allele/ID conflicts remain nullable. Valid sequence absent from the catalog is labeled \code{sequence_only}, with no invented ID. Generic opposite-event indels without a resolving known ID remain ambiguous.}
