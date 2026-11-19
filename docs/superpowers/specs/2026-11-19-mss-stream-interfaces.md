# Internal streaming interfaces

Existing exports and signatures stay intact. New private Rcpp exports:

    List mss_parser_batch_cpp(SEXP parser, RawVector input, int offset,
        bool eof, int max_rows, int batch_bytes, int field_count,
        IntegerVector keep)
    void mss_write_normalized_cpp(DataFrame table, SEXP connection,
        bool header, int block_bytes = 262144, bool native_encoding = false)

Reuse mss_parser_open_cpp(delimiter,max_record_bytes) with the new collector. Preserve mss_parser_feed_cpp's row-list contract. Collector mode must not silently misinterpret a parser already used through the other mode; reject incompatible mixing.

The new parser returns list(state, offset, nrow, columns). State is need_input when the chunk is exhausted and the partial native batch is retained; batch for a complete column batch including a final partial batch; done when EOF is drained. Offset is the zero-based next byte in the same raw chunk. R must neither slice/replay input nor double-write raw snapshot bytes. EOF uses an empty raw chunk and may be repeated until done. A nonempty EOF chunk is rejected. Non-batch states return nrow=0 and empty columns.

Header uses field_count=0, keep=integer(), max_rows=1. Its batch has one character vector of header names in columns[[1]] and nrow=1. Body field_count is the validated header width; keep is strictly increasing unique one-based source positions and may be empty. Columns is an unnamed list of character vectors in keep order, all length nrow. Row count is independent of selected-column count. Preserve decoded spelling until the existing single numeric conversion for mapped pos/effect/eaf. Validate and budget every source field, including unselected columns. Hold at most one complete overflow row and retain exactly historical batch membership.

The writer takes the existing R text-mode connection and validated normalized dataframe, emits the historical quoted header only when header=TRUE, emits every row, and returns invisible NULL. Bounded native buffers target block_bytes and stream oversized fields without constructing an unbounded encoded output row. Call base writeLines(block, connection, sep="", useBytes=TRUE) through the existing text connection only where native encoding matches the baseline. Keep the connection open for existing flush/hash/close. Same-platform bytes against the baseline sprintf/write.table path are the oracle. Retain the base emitter fallback for locale/encoding/platform cases whose exact parity cannot be established. Do not weaken snapshot, reference-mutation or checksum verification.

Round1 review ruling: the extra private native_encoding flag defaults to false, preserving the four-argument call as a conservative baseline-fallback call. Only the owned file producer may set it true, from the actual encoding captured when that same connection is created; the current global encoding option alone is not connection evidence. This is a revision of a new unshipped private interface, not a change to existing public MSS formals. Preserve the producer's historical output encoding by passing its captured value to file(), rather than forcing a new format. Do not inspect unstable R connection internals. Windows/locale/encoding guards still apply.

The native writer must decide whole-table fallback before emitting header or body. Marked inputs requiring whole-field translation (including Latin-1) use the allowed baseline formatter unless bounded conversion is implemented and proved. The native path may use the bytes of native/UTF-8 strings in a compatible UTF-8 locale without whole-field translation. Protect each CHARSXP, including an ALTREP-generated element or header name, for its entire use across allocating R callbacks. Regression tests must explicitly exercise the native flag on verified native connections, unknown explicitly re-encoding connections with the default flag, marked Latin-1 fallback, and a tiny forced-GC/deferred-character case. Do not let the original native-path tests silently become baseline-versus-baseline tests.

No lookup prefetch API is required. Optional integration needs a root-accepted concrete lookup contract; correctness always uses complete scalar demand.
