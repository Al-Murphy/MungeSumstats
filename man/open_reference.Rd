\name{open_reference}
\alias{open_reference}
\title{Open a verified reusable disk reference handle}
\usage{open_reference(reference, fasta = NULL)}
\arguments{\item{reference}{An \code{mss_reference}, a reference directory path, or an already open reference handle.}
\item{fasta}{Optional explicit runtime path to a byte-identical copy of the recorded FASTA. The original size/MD5 and all bundle artifacts, including FAI, are fully verified. No file is guessed or manifest rewritten.}}
\value{A process-local \code{mss_reference_handle} accepted by both harmonisation APIs.}
\description{Verifies reference MD5 checksums once and opens native file streams for reuse across batches and studies. Subsequent calls check file sizes and modification times.}
\details{The handle retains disk streams, contig metadata, and at most 28 MiB of reference/FASTA byte pages plus bounded page metadata, independently of total reference size. Keep the files immutable while it is live. Handles cannot be restored by serialization or reused in another R process; reopen the reference directory instead. Native streams are closed when the handle is garbage-collected.

For a transferred reference, use \code{open_reference("new/reference/path", fasta="new/reference.fa")}. The recorded original FASTA path need not exist. The returned handle's \code{fasta} is the actual normalized runtime path; its optional \code{fasta_binding} records both paths and the verified original size, FASTA MD5 and FAI MD5 under \code{mss-fasta-binding-1}. These physical binding details propagate to harmonisation provenance and streamed edition receipts only when explicitly used. Scientific checksums, source data, normalization versions and immutable reference-bundle bytes remain unchanged. Explicit rebinding always rereads the manifest and opens a fresh handle, including for a serialized or previously open handle input.}
\seealso{\code{\link{prepare_reference}}, \code{\link{harmonise_sumstats}}, \code{\link{harmonise_file}}}
