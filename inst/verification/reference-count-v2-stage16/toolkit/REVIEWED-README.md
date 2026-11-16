# Count-only native reference sizing

Operator diagnostic, pending independent review. No MSS API, default format,
normalizer, installed library or production reference is changed. The executable
only reads a plain FASTA and six-column prepared TSV and emits numeric JSON. It
does not call the writer or sorter, create output directories or authorize a build.

Build with a C++17 compiler against the exact MSS header recorded in PLAN.md:

```sh
c++ -std=c++17 -O2 -I /path/to/MungeSumstats/src count_reference.cpp -o count_reference
./count_reference /path/to/reference.fa /path/to/variants.tsv 3
```

An operator must authenticate and freeze input/core/executable pins, use the
resource guard, and recheck hashes/stat after counting. `qualify_prefix.py` is a
fixed, one-run qualification coordinator for the existing frozen incomplete1M
fixture; it is not a general full-input launcher. Input changes/errors refuse a
successful receipt. No public or untrusted-input authority is implied.

R is emitted record bytes; U/M raw sorter bytes/rows; K/N final key bytes/rows;
J rejected-file bytes; T the prepared input size copied later by the R wrapper;
C an explicit bound for FAI, version marker, contig aliases and operational
manifest. Incremental logical storage, after the prepared input already exists:

```
I = T + R + K + 8*N + J + C
P = max(I, R + 2*U + J + C, R + U + K + 8*N + J + C)
```

P accounts for raw inputs coexisting with merged output, and the last raw run
coexisting with final keys/offsets. It is a conservative logical envelope, not an
allocated-filesystem or memory certificate. Keep allocation/directory overhead,
concurrent reservations and the170GiB planning reserve separate. Recheck actual
free space and the160GiB hard floor; count both build orders before admission.

The counter uses the frozen native normalizer. Every record has a distinct byte
offset. Full sort-row equality includes that offset, so cross-record keys cannot
deduplicate; within-record distinct keys determine exact K/N without global sort.
Raw counts retain duplicate keys/aliases and unsupported/rejected occurrences.
Checked uint64 arithmetic refuses overflow rather than wrapping.

Differential tests independently inspect files from the real writer for formats2/3.
The standalone fixture writer is test-only. Ordinary and separate ASan/UBSan builds
cover29 cases plus overflow controls. Native file-byte parity is measured on
macOS; Windows diagnostic qualification is not claimed. Platform text-newline
accounting is explicit for the rejected text file.

The single existing1M prefix run reproduces all previously pinned native byte/row
counts in2.864s guarded count-process time. This excludes pre/post hashing and is
not an MSS whole-file harmonization speed result. The actual prepared production
inputs remain unsealed. There is no production count, full-index admission or
general genome/alias-fanout extrapolation.

## Coordinator branch extension

This branch adds `fai_upper_bytes`: exact FASTA contig-name bytes plus four
maximum-width uint64 numbers and native text framing per index row. This uses an
additional streaming header pass and preserves the production Fasta API. It is a
conservative small-file bound. Original R/U/K/N/J outputs remain exact and unchanged.
The original39-file diagnostic archive and2.864s prefix measurement belong to the
earlier frozen source; no new performance claim applies to this extension.

The sealed-input coordinator is implemented on this diagnostic branch and awaits
fresh independent review. See COORDINATOR.md for the trusted local profile, both-
producer seal gate, guarded complete worker and unconditional no-admission planning.
No full reference may be opened until both guarded production preparations finish,
independent review is GREEN, and root separately admits the pinned count-only profile.

The separate correction branch adds independently supplied expected publication
roots before control opens/mkdir, pre-read JSON size bounds and necessary native
response invariants. See COORDINATOR.md and the stage11 fixes report. These fixes
remain pending independent correction review; all no-admission gates remain.
