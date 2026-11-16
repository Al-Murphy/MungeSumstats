# Full-reference capacity planning for reviewed typed format 3

**Decision: no full-index admission yet.** Reviewed opt-in format 3 improves the
storage shape, but both complete history-expanded inputs remain unpublished and
their final capacity is unknown. Neither “both fit” nor “the first cannot fit” is
established. Keep format 2 as default, the existing 8GiB memory/160GiB free-space
guards, and a **170GiB planning reserve**. This assessment launches no full index
and performs no cloud, cleanup, source-preparation or production-library operation.

## Evidence and current state

The storage implementation is frozen at
`ece4fcdfad4aa9c2ce85434f67694190872442dd`, core SHA-256
`9c06526ec9185940d36e3002ff5d9c06a395d34e6a29d1cddb4b7ea119480011`.
It uses reviewed A7/S5 records and unchanged scientific implementation 7. See
`mss-context-format3-independent-review-stage7.md` for that separate review.

At **2026-11-16 09:14:17 UTC**, a metadata-only snapshot reported
**479,232,000,000 bytes free**. Subtracting 170GiB = 182,536,110,080 bytes leaves
**296,695,889,920 bytes**, before further preparation growth, other work or an
allocation allowance. The active guards still reported `running`; neither final
published `variants.tsv` existed. Growing intermediate `sites.tsv` sizes were
27,869,754,883 bytes for GRCh38 and 17,523,082,865 bytes for GRCh37. These are
occupied bytes, not final-input predictions or lower bounds on final TSV sizes:
identifier disposition and complete historical expansion are unfinished.

Only small live guard JSON and filesystem metadata were read from active work.
No growing TSV/SQLite/VCF was scanned, hashed, paused, moved or modified. The new
streaming analysis read the eight existing frozen 100k windows and the frozen 1M
prefix: 56,703,923 source bytes and 1.8M row occurrences, with overlapping scope
between the prefix/windows. It created no indexes. The unchanged operator guard
completed successfully in 13.31s, with **21,741,568 bytes sampled peak memory** and
480,248,918,016 bytes minimum free space. The first sandboxed attempt refused
closed because `ps` was unavailable; its censored receipt is retained. The retry
used process-metadata permission and the same 8GiB/160GiB bounds.

## Exact retained and temporary byte accounting

Let:

- `T`: completed prepared input TSV bytes, copied into the new index.
- `R`: all emitted typed A and S record bytes.
- `U, M`: bytes and row count of raw sort key rows, before any deduplication.
- `K, N`: final sorted/deduplicated key bytes and row count.
- `J`: rejected-row output bytes.
- `C`: FAI, 19-byte format marker, explicit contig-alias table and operational
  manifest allowance. The external shared FASTA is accounted separately.

Relative to free space measured after the prepared input already exists:

```
I = T + R + K + 8*N + J + C
P <= max(I, R + 2*U + J + C, R + U + K + 8*N + J + C)
```

`I` is retained index size. `P` is a conservative logical phase envelope. During
multilevel 32-way merges, input runs coexist with the new merged output. All live
inputs are bounded by total raw bytes `U`, and new output adds at most another
`U`. During final key/offset emission the last raw run remains until the read
finishes; final keys and the 8-byte offset table accumulate beside it. Counting
only final keys, or assuming sort inputs disappear while output is written,
understates scratch. Source copy and manifest publication happen after native
sorting; the retained phase accounts for that later `T` copy.

These are logical bytes, not a universal allocated-filesystem ceiling. Admission
must additionally reserve an explicit allocation/directory/metadata allowance,
snapshots or other reserved work, and the planning reserve. Actual free-space
measurements and the hard guard remain necessary. Do not infer available space
from the sum of files proposed for deletion or APFS shared extents.

## Conservative source-only upper bound

The new bounded helper reads **every source row**, including unresolved/symbolic
contexts, unsupported ALT strings, nonpositive/invalid positions and all primary
IDs/aliases. It never filters or silently skips a row to obtain a smaller bound.
It verifies source SHA-256 and unchanged stat before accepting a result.

For row `i`, define byte lengths `c_i=CHR`, `r_i=original REF`, `L_i=complete ALT
field`, `b_i=primary ID`; `n_i` is its raw comma-split ALT token count, including
repeated/trailing empty tokens. As in the C++ splitter, a completely empty ALT
field has zero tokens. Let `p_i` indicate a nonmissing primary ID, and let its
nonmissing alias tokens have byte lengths `a_ij`. Missing means exactly empty,
`.` or `NA`, matching the writer. Repeated aliases still count in raw sort rows.

An S record has four tabs, its one-byte tag and LF: six framing bytes. Its signed
original position occupies at most 20 bytes. An A record has eight framing bytes
and two positive signed-64-bit positions, at most 19 bytes each. Trimming and
left reanchoring cannot increase canonical REF or ALT length beyond their own
original allele lengths. Therefore:

```
Rmax = sum_i (26 + c_i + r_i + L_i)
     + sum_i sum_ALT (46 + c_i + r_i + len(ALT) + b_i)

Umax = sum_i (c_i + 45)
     + sum_i n_i * [2*(c_i + 44)
                    + p_i*(b_i + 24)
                    + sum_j(a_ij + 24)]

Mmax = Q + sum_i n_i * [2 + p_i + alias_token_count_i]
Jmax = 11 + 54 * sum_i max(1, n_i)
```

`Q` is source-row count. Every possible S is included, even an invalid-position
row which the writer would reject without S. Every ALT token is treated as a
possible A, including unsupported tokens which actually produce no A. This is an
upper bound, not an asserted resolved-variant count. The rejection bound covers
one invalid-position rejection even with no ALT tokens, or one rejection per ALT;
54 bytes allow a 20-digit row number and a status longer than those currently used.

Each sort-row bound includes a **20-digit unsigned record pointer**, tab and LF.
S coordinates allow 20 signed digits; C coordinates allow 19 positive digits.
There are two raw C keys per possible A, even when original/canonical positions
coincide and final emission later deduplicates them. Alias/current-ID I rows are
counted for every potential A occurrence. Byte lengths, rather than character
counts, preserve Unicode alias accounting. The model also checks the unsigned
64-bit address envelope; exceeding it is a refusal, not a wraparound assumption.

Since `K <= Umax` and `N <= Mmax`:

```
Iupper = T + Rmax + Umax + 8*Mmax + Jmax + C
Pupper = max(Iupper, Rmax + 2*Umax + 8*Mmax + Jmax + C)
```

The supplied model reserves 1MiB for the small operational RDS metadata, plus
explicit FAI/marker/alias-table bytes. Final deployment metadata must be checked
against that allowance; additional artifacts increase `C`. Unbounded future
metadata, disk allocation or concurrent writers are not assumed away.

The bound deliberately trades precision for a cheap pass without normalization.
For the 1M prefix it is 447,501,755 bytes retained and 618,255,719 bytes phase
envelope, versus approximately 198MB/268MB exact native accounting. A conservative
denial would not prove physical impossibility. Use the exact count-only option
below before changing guards, reducing content or provisioning from a ratio.

## Calibration: actual prefix and frozen windows

No new sample index was built. The helper re-expresses every immutable format-2
canonical record as A7/S5, tracking its new physical byte offset and all raw keys.
Source/current-ID/original-REF/full-ALT pairing is checked at every occurrence.
The 1M recount reproduces the existing format-3 record/key/offset file sizes
exactly. Other windows are **exact counterfactual byte counts**, not measured
format-3 builds. Operational manifests are excluded from the native column.

| Frozen input | Source bytes | Typed-3 native retained | Modeled typed-3 phase envelope |
|---|---:|---:|---:|
| GRCh38 0% | 2,998,568 | 18,487,414 | 24,594,379 |
| GRCh38 30% | 3,079,940 | 18,986,184 | 25,298,486 |
| GRCh38 60% | 3,259,335 | 19,857,052 | 26,413,782 |
| GRCh38 90% | 3,197,356 | 19,458,518 | 25,841,643 |
| GRCh37 0% | 2,994,390 | 18,393,450 | 24,457,308 |
| GRCh37 30% | 3,300,549 | 20,215,505 | 26,903,179 |
| GRCh37 60% | 3,269,195 | 19,805,500 | 26,310,306 |
| GRCh37 90% | 3,196,206 | 19,912,892 | 26,611,971 |
| Actual 1M prefix | 31,408,384 | 198,422,813 | 267,719,951 |

The prefix has 1M S occurrences, **1,242,071 valid A occurrences**, and 1,242,077
raw ALT tokens; the upper model retains all six non-emitted tokens in its bound.
It has 4,726,213 raw sort rows versus 3,484,369 final rows, with exact typed-3
`R=63,883,154`, `U=100,705,522`, `K=75,230,231`. Adding its actual 743-byte
manifest yields **198,423,556 public retained bytes**, versus 232,818,833 for 2:
the reviewed **14.7734%** prefix reduction.

The phase envelope is **267,719,951 bytes**, larger than the author's sampled
257,008,054-byte work-directory peak. Samples can miss the end of final emission;
the smaller sample is not an admission ceiling. Old-2 modeled accounting was
302,403,057 bytes; typed-3 modeled reduction here is 11.47%, distinct from retained
reduction. The prefix measurement establishes no RAM improvement.

The eight windows have typed retained/source ratios 6.058–6.230 and modeled phase
ratios 8.048–8.326. Their retained savings versus 2 range 13.445–15.259%. **All nine
inputs have zero historical alias tokens**. They cannot predict the unseen
remainder or final history-expanded inputs. The model bounds all nine recounts,
plus existing frozen 256/2,048/10,000-ALT and signed/invalid-position/duplicate/
Unicode fixtures. Small empty/signed-field arithmetic controls also pass.

## Why source bytes or a universal percentage are insufficient

Typed 3 removes repeated full ALT lists from A records, but still stores canonical
REF and current ID **once per valid ALT occurrence**. The source-only bound thus
retains the `sum(n_i*r_i)` term. An allocation-free shape calculation with a
1,000,000-base REF and 10,000 distinct eight-base ALTs has about 1.09MB of source
row text yet over **10GB of possible typed record bytes**. Such untrimmed events
are possible when REF and ALT have differing first/last bases. No such huge
fixture was allocated or built. This proves why source-GB multiplication and the
14.77% prefix saving cannot be universal capacity certificates.

Historical expansion is also representation dependent. Adding one alias to an
existing row changes no A/S payload, but adds one raw I key per ALT occurrence.
The current preparer instead emits **another complete context row for each
verified historical alias**. Each repeated row adds another S, all valid A
occurrences, original/canonical C keys, primary I keys and the historical I keys.
Rows remain distinct physical occurrences even when their canonical events agree.
Global `verified_alias_pairs` alone is insufficient: the same identifier can
participate in multiple site contexts with different ALT fanout. Count the final
emitted TSV with its actual ID/alias disposition. The reviewed preparer is pinned
at `84199d3b37273ab56c838aac385ffcdea40197bc`; this assessment changes none of its
existing identifier policies and adds no filtering or quarantine exceptions.

## Exact count-only refinement without full index scratch

After both inputs are sealed, a separately reviewed dry-run writer can reuse the
frozen normalizer and exact writer predicates, emit no index files, and count
`R,U,K,N,J` plus maxima. This can determine exact retained logical size without
building or globally sorting key rows:

1. Count each exact typed record's length and advance its virtual byte offset.
2. Enumerate all of that record's raw C/I/S keys with that decimal pointer for `U`.
3. Deduplicate keys **within the record** for exact `K,N`.

The third step is sufficient because equality of full sort rows includes the
physical pointer. Distinct nonempty emitted records have distinct offsets;
therefore keys from different records cannot deduplicate. Sorting changes their
order, not byte sums. This proof includes repeated historical rows, duplicate
aliases and duplicate original/canonical C keys. The new recount exercises that
identity and reproduces the actual prefix files. A production normalize/count-only
CLI is **not implemented by this report**; existing frozen canonical records are
used only for calibration. Its future implementation must be tested against exact
writer payloads, retain all raw counts and fail closed on errors or mutation.

## Next admission procedure

Finish both guarded preparations first. Require immutable input digests/sizes,
complete history/source receipts and final emitted-row counts. Run a bounded
source-count pass on the **final** TSVs, record REF×ALT/ID/alias-weighted sums,
maximum line/key/REF/ALT sizes, and compare `Pupper` with actual free space. The
helper defaults to 64MiB per input and 8MiB per line for these frozen experiments;
these are operational analysis limits, not permitted variant filters. A larger
real line/input must trigger a reviewed limit increase or a more streaming parser,
never truncation or row omission. Use exact normalization/count-only refinement
if the conservative envelope is too loose.

For order A then B, let `I_A,I_B` include every new retained index artifact and
`P_A,P_B` the guarded phase envelopes. Relative to free space measured **after
both prepared inputs and needed FASTAs already exist**, the requirement is:

```
max(P_A, I_A + P_B, I_A + I_B)
  + allocation/metadata allowance
  + separately reserved concurrent work
  + 170GiB planning reserve
```

Evaluate both build orders and choose the smaller requirement. Admit one index
at a time, verify it completely, then remeasure actual free space before B. Keep
the 160GiB hard floor and 8GiB process guard throughout construction, hashing and
verification. Disk sizing does not prove resident-memory compliance: maximum
source/key lengths and chosen sort-buffer rows matter, and ordinary R/package
startup memory remains. A smaller sort buffer can reduce memory without dropping
content; the same `2*U` disk envelope still applies.

Both runtime closures are needed for complete cross-build projection: all
manifest-listed files, exact `source.tsv`, operational manifest, and the external
plain seekable FASTAs (currently 3,339,671,278 and 3,275,274,505 bytes). Exclude
already-present FASTAs from incremental index bytes, but include them when sizing
a new workspace. Object storage alone is not the supported seekable runtime.

If capacity is insufficient, consider only verified offload of **sealed inactive**
artifacts with recoverable receipts and end-to-end content verification, then
measure actually reclaimed free space. Do not remove active acquisitions,
intermediate sites or SQLite state. If both closures/build scratch still do not
fit, size a separately authorized workspace from the sealed counts, preferably
US East/US Central/Canada East, and obtain a current budget decision before any
resource operation. No provider operation or new funding assumption is made here.

## Reproduction and pins

New evidence is under `results/reference-format3-capacity-stage7/`:
`model.py`, `inputs.json`, `model-v2.json`, `summary.json`, `check_model.py`,
`model-checks.json`, `snapshot.json`, both guard receipts/logs and `pin-index.json`.
The input manifest pins all nine frozen source byte streams and their calibration
indexes. `pin-index.json` pins model/guard/reference/preparer source and outputs.

The successful invocation is recorded verbatim in `guard-v2.json`. Reproduction
uses a fresh output path; it must not overwrite receipts. The helper rejects
changing source metadata/digests, over-limit lines/files and `-pending-` source
paths. Every result explicitly says `full_reference_admitted=false` and
`index_builds_executed=0`. Full cardinality, final historical growth and final
capacity remain unknown until the completed inputs are sealed and counted.
