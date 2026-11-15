\name{open_reference}
\alias{open_reference}
\title{Open a verified reusable disk reference handle}
\usage{open_reference(reference)}
\arguments{\item{reference}{An \code{mss_reference}, a reference directory path, or an already open reference handle.}}
\value{A process-local \code{mss_reference_handle} accepted by both harmonisation APIs.}
\description{Verifies reference MD5 checksums once and opens native file streams for reuse across batches and studies. Subsequent calls check file sizes and modification times.}
\details{The handle retains disk streams and contig metadata, not reference variants or FASTA sequences. Keep the files immutable while it is live. Handles cannot be restored by serialization or reused in another R process; reopen the reference directory instead. Native streams are closed when the handle is garbage-collected.}
\seealso{\code{\link{prepare_reference}}, \code{\link{harmonise_sumstats}}, \code{\link{harmonise_file}}}
