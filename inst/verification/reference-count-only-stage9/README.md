# Reviewed count-only reference capacity diagnostic

Independent review is **GREEN** for the isolated operator diagnostic at
`1bf1a08685488798496e9d7978b1792dbd290a9d`. This adds no MSS API, default,
normalizer, reference format or installed-library change. It grants no production
input/index admission. The [source and build procedure](source/README.md) and
[independent review](independent/reference-count-only-independent-review-stage9.md)
are included.

The C++17 executable reads the frozen native normalizer's inputs and counts exact
record bytes, raw sort bytes/rows, final key bytes/rows and rejected bytes. It
creates no reference index, sort runs or FAI. Each complete lookup key includes
its physical record offset, so within-record deduplication yields exact final
byte/count totals without global sorting. Repeated historical rows remain distinct;
duplicate aliases/coordinate keys are retained in raw sort totals.

Fresh review verifies all10 source pins and the unchanged production header,
reruns29 author cases and adds22 actual-writer parity/refusal cases. Overflow
controls pass. Separate strict ASan/UBSan evidence covers56 successful native
invocations with empty stderr plus overflow control. The malformed diagnostic
version parser RED/fix and earlier sanitizer-runner inspection limitation remain
preserved. No Windows diagnostic sizing claim is made; platform text-newline
handling is explicit.

One run on the existing incomplete1M-site prefix matches every previously pinned
R/U/K/N/J value exactly. Count-process time is2.864s with6,389,760B sampled peak RSS;
this excludes streaming pre/post hashing and is **not** an MSS whole-file speed
benchmark. No index was rebuilt. This prefix has no historical alias tokens and
cannot predict complete history-expanded production inputs.

The retained/temporary envelope includes source copy, coexisting merge inputs and
outputs, and explicit small-file bounds. Filesystem allocation, concurrent work,
both build orders and170GiB planning reserve remain separate. See
[candidate report](reports/reference-count-only-stage9.md) for the measured scope.
The fixed prefix coordinator is not a full-input launcher. A future full count
requires sealed authenticated input closures, explicit allowances/guards and fresh
before/after pins; active GRCh37/38 preparations were untouched.

`archive-map.json` binds original/archived sizes and SHA256; logs use lossless gzip.
Original indexes remain verbatim. `omissions.json` pins generated binaries and
external inputs. All artifacts were scanned against actual private test/account
values; credentials, account inventories and hook memory are absent. Every indexed
artifact is explicitly checked as tracked before transfer.
