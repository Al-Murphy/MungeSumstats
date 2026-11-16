# Independent review: opt-in typed reference format 3

**Verdict: GREEN for the documented opt-in typed reference-format contract of frozen `ece4fcdfad4aa9c2ce85434f67694190872442dd`.**

No material correctness finding remains in this reviewed candidate. This approves
the frozen opt-in storage implementation for subsequent installed-package and
native operating-system CI gates. It does not approve making format 3 the default,
admit either complete production reference to a capacity budget, or claim that
full installed-package/native Windows CI for format 3 has already passed.

## Frozen scope and identities

The reviewed clone is `/private/tmp/gwas-harvest-stage2/mss-context-format3`, branch
`perf/typed-reference-format3`, scientific base
`5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`. Its production worktree was clean on
the final check. All independent compiles use a fresh `git archive` of the frozen
commit under `review-mss-context-format3-stage7/source`; no production source,
library, reference preparation or cloud resource was modified.

| Artifact | SHA-256 |
|---|---|
| `src/preservation_core.h` | `9c06526ec9185940d36e3002ff5d9c06a395d34e6a29d1cddb4b7ea119480011` |
| `src/preservation.cpp` | `6d7281661962b2db9f74ae53569841336d810315ab417878b51969d82f63a49e` |
| `src/RcppExports.cpp` | `2b0a3d1bbf1809120dd60d1b1c802b4355f958f2547e545262909feb0db5a1d3` |
| `R/RcppExports.R` | `895594e02ee6c638600abe6085e51d1b0243e706a8dbdccd0d353e9de3b568f2` |
| `R/preservation.R` | `60e1936f2e5378b534ed387351412cb46eff01b58ed6f18efd950f2f83a4dde1` |
| `R/projection.R` | `920962a89afab78315ed48506c4fae5f4ae33e5bbdc4d79d173838a770b89531` |
| `inst/PRESERVATION_CONTRACT.md` | `c2a9eee134e567ebd5e5bcdd5dc02e7ddeacf67cb347bae8ba13ff451c9f56e8` |
| Installed prototype DLL | `d4f43e62b77a2135f94084ac3fe15105cbc6809fed388d2f4e7bc4ee9bc926da` |

The installed DLL used by the small public R checks is under
`results/mss-context-format3-stage7/installed-lib/MungeSumstats/libs/`.
The complete C++ suffix beginning at `Result harmonise(` is byte-identical to
scientific implementation 7, SHA-256
`47c1a1261aaff1421293d41089df74fedba37f264ed71ae78e930a9030b90d47`.
Reference storage version 3 is distinct from unchanged normalization 7,
delimited parser 1 and projection 2.

## Preservation argument and limits

The A7/S5 split preserves the information needed by scientific 7. Every parsed
numeric-position source occurrence has a separate S physical pointer; every valid
ALT occurrence has a separate A pointer. Duplicate source rows and repeated ALT
occurrences remain separate for raw lookup caps. Repeated alias keys pointing to
the same physical A are deduplicated as before. C/I raw key multiplicities are
unchanged; raw source bytes retain all unsupported records and identifiers.

An A record contains canonical event, current primary ID and original position.
Current-ID closure, event closure, main-pool intersection and candidate count use
these A-derived fields. S IDs therefore carry no omitted scientific evidence.
Unsupported-only contexts remain S records and remain discoverable through their
site position; neither format invents an I mapping for an unsupported-only row.

Each S reconstructs its own chromosome, original position, original REF and empty
canonical ALT before semantic merging. Consequently REF spans cannot be pooled
across unrelated contexts. Full ALT payload is loaded once through S lookup;
C/I processing does not repeatedly materialize it. Same-span ALT union remains
the existing scientific-7 operation. A singleton S keeps its original raw token
occurrences/order, including duplicates charged to the work budget.

Format 3 does not broaden scientific 7's eligible universe into every overlapping
unlinked genomic record. Its current-ID/exact-event closure, missing-position
complex-allele withholding, composed-background necessary-condition deferrals and
trusted declared REF/ALT contract are unchanged. `source_allele_context_required`
is a conservative metadata-required result, not proof that two actual mappings
exist. The prior scientific report and older RED reports retain their own scopes.

## Independent verification actually executed

| Gate | Result |
|---|---|
| Original scientific R suite, format 2 | 65 blocks / 1,170 assertions |
| Original scientific R suite, opt-in 3 | 65 blocks / 1,170 assertions |
| New format API suite through independently compiled source | 7 blocks / 65 assertions |
| Same API suite through pinned installed DLL | 7 blocks / 65 assertions |
| Additional installed public API, receipt and integrity controls | 41 passed |
| Finite SNP mathematical oracle, each format | 16,128 rows, zero violations |
| Haplotype/contract mathematical oracle, each format | 79,904 rows, zero violations |
| Explicit scientific/resource controls, each format | 297 rows: 156 resolved, 132 nullable, 9 refusals |
| Direct complete native output equality, 2 versus 3 | 96,032 mathematical-oracle rows byte-exact |
| Old-7 versus new default-2 writer matrix | 6 source orders × 6 buffers; all six native payloads byte-exact |
| Typed-3 writer matrix | Same 36 combinations; exact raw key multiplicities |
| Independent selected-record structural probes | 73 cases: 18 accepted positives, 55 refused corruptions |
| Extra shared-event/current-ID ordering controls | 6 source orders × 4 queries × 2 formats; exact expected fields |
| Explicit work-order qualification | 8 format/order/cap combinations matched documented behavior |

The two R science runs use unchanged scientific expectations. For opt-in 3, the
prepare default is injected in the review process and exactly one old provenance
literal is changed from format 2 to 3 in a private test copy. Source tests and the
author clone remain untouched; the adaptation has its own receipt.

The SNP oracle covers minimal, left-, right- and both-padded representations,
allele orders and strand/effect orientations. Its 4,608 resolved payloads per
format and 11,520 uncertainty payloads match independent expectations. The larger
haplotype corpus checks 2,104 independently valid resolved payloads, 17,736
independent ALT–ALT explanations, unknown-anchor withholding, composed background
guards and coordinate/ID conflicts. The 297 explicit controls include both earlier
minimal background counterexamples, repeat reanchoring, ID/event cycles, foreign
positions, unknown markers, duplicate ALT tokens, signature/normalization work
budgets and positive biallelic indel controls. These are mathematical expectations,
not old-output gold.

The additional ID-order controls check that collapsing several current IDs for
one event leaves the returned canonical ID empty, while direct unambiguous ID
queries retain the correct ID. Changing decimal physical offsets did not leak a
first-encounter primary ID. No new canonical resolution or effect-direction gain
is claimed.

## Strict reader and public integrity boundaries

Independent malformed-index probes bypass R checksums to exercise selected native
records: tags/field width, canonical DNA/distinct alleles, UTF-8 overlong/surrogate/
out-of-range/truncated encodings, NUL, empty fields, canonical integer spelling,
signed/unsigned overflow, wrong C/S coordinates/types, missing LF, offsets into
records, missing/extra key fields, truncated binary offset tables and malformed
markers. Positive cases include opaque Unicode primary IDs/aliases, canonically
different Unicode aliases, long records/keys and UTF-8 crossing 64KiB boundaries,
both original/canonical C anchors, symbolic/empty ALT contexts, and signed S
positions including the minimum signed 64-bit integer.

The duplicate-occurrence fixture independently verifies S and A physical counts,
raw I/S cap refusal boundaries and exact raw key multiplicity across formats.
Singleton raw-token work refusals remain equal in both formats; the early
prototype's raw-token deduplication is absent from this frozen code.

Selected-record validation is not global certification of unread native records,
key sort order or alias membership. An A row cannot alone prove that an arbitrary
I alias belongs to its primary ID. Public reference open retains immutable whole-
artifact checksum verification, while source bytes preserve the alias relation.
Our public R tests independently refused a corrupted unqueried payload, mismatched
marker/manifest, wrong layout/order metadata and invalid UTF-8 source publication.
No partial invalid-source output was published. Arbitrary unverified native files
outside this selected-record/checksum contract are not certified.

## Compatibility, receipts and resource order

Default reference output remains 2. All six native payloads were byte-identical to
frozen 7 in the 36-case writer matrix, including a 300-alias row forcing sorter
merge levels. The installed native R wrapper still accepts four supplied arguments
and defaults its fifth to 2; registered arity is 5. This preserves wrapper source
compatibility, not an undocumented direct four-argument `.Call` ABI. The old native
reader explicitly refuses format 3.

All four source/target projection-format combinations produce equal scientific
rows and identify their formats honestly. Complete 2→2 projection receipts are
identical to the previous scientific-7 R function. Combinations containing 3 add
source/target format fields; unlike formats use `mixed`. Public normalized/source
rows, aliases, contig aliases and version receipts also pass installed checks.

The declared context-order change is real and reproduced independently:

| Context order / cap | Format 2 | Format 3 |
|---|---|---|
| Small competing context first / 5 | ambiguous | ambiguous |
| Long noncompeting context first / 5 | work refusal | ambiguous |
| Either order / 200 | ambiguous | ambiguous |

Format 3 sorts merged S contexts by chromosome/original position/original REF;
these fields uniquely identify a merged context. Raw record caps remain exact.
This is a versioned resource-order difference, not unconditional exception/status
parity. Every successful result in this fixture withholds canonical identity,
effect, EAF and orientation; there is no resource-exhaustion canonical promotion.

## Existing bounded measurements independently corroborated

This reviewer did not rerun performance measurements. Existing frozen harness,
DLL and memory-helper pins were checked. All 15 author heavy-command guard
receipts report success within 8GiB sampled tree/physical bounds, 160GiB free-space
floor and 4GiB aggregate experiment storage cap. Existing small output artifacts
and actual retained index file sizes were rechecked independently.

The actual first-million-site indexes total 232,818,833 bytes for 2 and
198,423,556 bytes for 3: **14.7734% retained reduction**, including source copy and
manifest and excluding the shared external FASTA. This is not complete-reference
capacity admission. No RAM improvement is established by this prefix.

The 256/2,048/10,000-ALT author artifacts retain all source occurrences, have exact
19-field public normalized/source equality, and withhold every canonical payload
field with `source_allele_context_required`. Their actual file sizes match the
reported totals. These are complete conservative generic-query paths with stated
proof caps, rather than resolved-variant throughput. The high-ALT gain is credible
for the repeated-payload shape removed by A/S factoring; it cannot be generalized
to ordinary biallelic GWAS or corpus throughput. Runs are single fixed-order pairs
on a shared host with uncontrolled OS cache. Disk and process-tree peaks are
sampled; they are not exact filesystem-event maxima. Small timing differences do
not establish general CPU gains.

## Evidence and remaining release gates

Independent evidence root:
`/private/tmp/gwas-harvest-stage2/review-mss-context-format3-stage7/`.
Primary receipts: `source-pins.json`, `science-adaptation.json`, `science.log`,
`structural-summary.json`, `compatibility-summary.json`, `installed-summary.txt`,
`installed-controls.tsv`, `final-controls-summary.json`, the per-format oracle
summaries/outputs, and `evidence-index.json`. Private review-harness compile/path/
JSON-shape mistakes are retained in earlier review logs; corrected scripts and
final successful logs supply the gates above. They were not production failures.

Remaining coordinator gates are the full installed package suite and actual
native multi-OS CI for this exact format-3 candidate, followed by a separate
complete-reference capacity decision before any production/default change. This
review neither modifies the active reference preparations nor weakens their
resource guards.
