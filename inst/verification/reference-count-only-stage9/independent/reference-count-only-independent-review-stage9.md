# Reference count-only diagnostic — independent review

**Decision: GREEN for the isolated count-only diagnostic; no blocking finding.**
This admits neither a complete input nor a production reference build. No prefix
1M/full-input count, active preparation read, index qualification or timing rerun
was performed by this reviewer.

Frozen source: `1bf1a08685488798496e9d7978b1792dbd290a9d`, tree
`7cce0bd3043555319835c3d878b631ab8013f30a`, independently cloned to
`/private/tmp/gwas-harvest-stage2/reference-count-independent-stage9`.
All ten tracked source files match the author index, and the reviewer clone
remains clean. Author evidence index SHA256
`94b45796456938b95509fcb47cbf7a0fb6ec01ec626e5b71c0b81b7d6fdf84b4`
and every listed artifact/generated-binary hash and size were verified. The
included production header remains exactly
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`;
it was read only, not changed or copied over production source.

## Exactness proof and writer parity

The proof is sound for the frozen writer. Each emitted S/A record has positive
byte length and therefore a distinct record offset. A raw sorter row is exactly
`key + TAB + decimal_record_offset + LF`; full-row equality includes the offset.
Different physical records cannot deduplicate even if all scientific fields and
keys repeat. Consequently the sum of distinct keys **within each physical record**
is the final K/N after the writer's global full-row deduplication. The counter
preserves raw duplicate coordinate keys and primary/historical identifiers in U/M.
It never deduplicates repeated historical source contexts across records.

The counter uses the same frozen FASTA/normalizer and source parsing/disposition
as `prepare()`: strict six-field/header handling, version3 UTF8/framing validation,
`stoll` parse consumption, uppercase allele tokens, independent site contexts,
normalized A records, rejected occurrences and primary/alias missing-value rules.
Format2 remains the writer's permissive byte-text behavior; version3 rejects
invalid UTF8/NUL. Counter position/byte formatting matches the emitted native
record layout, including decimal offsets crossing digit widths. The counter adds
checked uint64 arithmetic and refuses overflow rather than emitting wrapped counts.
Its source read error check is fail closed, and successful numeric JSON is emitted
only after traversal. It does not invoke prepare/sorter/output-directory creation.

Fresh local C++17 builds succeeded and produced the same binary hashes as the
author's three ordinary binaries. The unchanged **29 differential cases passed**,
and the actual addition/stream overflow control passed with empty stderr.
An additional **22 small version2/3 cases** compare actual writer file sizes and
refusal behavior. They cover repeated normalized ALT events, shared keys across
111 records and offset widths, int64 position endpoints, CRLF FASTA, interior-CR
aliases, NUL/invalid/overlong/surrogate UTF8, truncated/extra TSV fields and empty
input. Successful cases match records/keys/offsets/rejected bytes; invalid cases
refuse on both sides without counter JSON. Fixture input hashes and directory
contents are unchanged by the counter, with no counter-created index.

These tests validate byte counting, not scientific reference completeness or
readability of every permissively accepted legacy format2 input. The raw ALT-token
diagnostic counts source tokens even for an invalid-position row; emitted A/S and
rejection counters still follow their explicit dispositions. No token count is
silently substituted for successful allele records.

## Phase envelope and platform limits

The phase formula includes retained source copy T, R, K, 8*N, J and explicit C.
The frozen R wrapper creates `source.tsv` **after** native prepare/sort returns,
so T need not coexist with raw sort runs; it is included in retained I. During
merge, input runs plus the currently written merge output are bounded by 2U.
During final emission the last raw run can coexist with K/offsets. Applying the
stated C allowance conservatively to every phase is appropriate. FAI, marker,
contig-alias and manifest bounds remain coordinator obligations; filesystem
allocation/directory overhead, external existing inputs, concurrent reservations,
both build orders and free-space reserves remain separate.

Binary records/keys use LF. The rejected file is text-mode in the real writer;
the counter explicitly adds Windows newline expansion only for that stream.
Native Windows newline/filesystem behavior has not been independently executed
here and no Windows sizing qualification is granted. Memory is bounded by the
core page cache plus current line/alleles/aliases and per-record unique-key state,
not a fixed global RAM certificate for arbitrary single-row fanout.

The author's strict ASan/UBSan evidence (56 successful native calls with empty
stderr plus overflow control) was hash-verified, not rerun by this reviewer.
The already measured incomplete1M timing excludes pre/post hashes and is not an
MSS benchmark. No extrapolation, input-sealing approval, full coordinator or unsafe
admission boolean has been added. The frozen prefix coordinator remains explicitly
prefix-specific; a future full-input coordinator needs authenticated closure,
before/after pins, guards and fresh explicit sizing allowances.

## Evidence

`/private/tmp/gwas-harvest-stage2/results/reference-count-independent-stage9/`
contains `run_review.py`, `verification.json`, three fresh build logs,
`author-29-fresh.log` and `independent-small-cases.json`. The companion
`evidence-index.json` binds these files and this report. Logs/evidence were created
under umask077. No credentials, private account environment, provider/network
operation, source edit, main merge or production/default change was used.
