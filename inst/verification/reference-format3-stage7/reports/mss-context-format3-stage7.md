# MSS opt-in typed reference format 3 prototype

The isolated prototype is frozen at **ece4fcdfad4aa9c2ce85434f67694190872442dd**
in `mss-context-format3`, based on scientific implementation 7 commit
`5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`. **Default reference format remains 2.**
Small correctness gates and the authorized bounded measurements are complete.
The compact representation reduces actual-prefix storage and removes repeated
full-site payload from the synthetic high-ALT stress cases. Complete-reference
capacity admission and default adoption remain unapproved.

## Public interface and physical layout

The appended R option is
`prepare_reference(..., reference_format=c("mss-preservation-2", "mss-preservation-3"))`.
The native R wrapper appends `format_version=2L`; its registered `.Call` arity is
five. Existing four-argument R wrapper calls continue using format 2. The portable
C++ writer similarly appends `int format_version=2`. A locally installed package
was used to verify registration and the public API, independently of `sourceCpp`.

Format 3 retains the same immutable artifact filenames, with these typed records:

```
A<TAB>chr<TAB>canonical_pos<TAB>REF<TAB>ALT<TAB>current_id<TAB>original_pos
S<TAB>chr<TAB>original_pos<TAB>original_REF<TAB>complete_ALT_list
```

Every parsed numeric-position source row occurrence gets a distinct S offset;
every successfully normalized ALT occurrence gets a distinct A offset. Identical
raw occurrences remain distinct. Sort deduplication removes only repeated keys
pointing to the same offset. Unsupported/nonsequence ALT contexts remain in S,
including rows with no valid A. Unparseable positions keep the existing
`rejected.tsv` policy and exact original `source.tsv`, without inventing an S
position. Signed nonpositive S positions preserve existing format-2 behavior;
canonical A positions must be positive.

A holds the current primary identifier; aliases remain exact opaque I keys.
C/I lookup never loads or unions the full site ALT list. Original positions locate
S. S reconstructs chromosome, original position, original REF and empty candidate
ALT **before** same-span merging, preserving independently paired REF spans.
Scientific current-ID closure remains A-derived: omitted S IDs do not participate
in canonical identity selection, closure, ALT–ALT interpretation or EAF checks.
All original identifiers, aliases and formatting remain byte-exact in `source.tsv`.
There is no first-A-as-S reuse.

Readers accept 2 and 3; old format-2-only readers reject 3. The data marker and R
manifest identify 3 separately from unchanged scientific implementation 7,
delimited parser 1 and projection algorithm 2. Any projection combination containing
3 records explicit source and target reference formats; common `reference_format`
is `mixed` when they differ. All four combinations are tested.

## Structural integrity and resource order

Format 3 validates complete selected A/S records and binary-search key rows:
field counts, tags and compatible key types, coordinate agreement, concrete and
distinct canonical A DNA, canonical signed/unsigned integer spelling and range,
UTF-8 without NUL/TAB/LF framing, in-range physical line-boundary offsets and a
terminating LF, including cache-page-spanning records. UTF-8 identifiers are not
case-folded, normalized or silently discarded. Alias membership is not derivable
from A alone; checksummed keys and the complete source retain that relation.
Unqueried unrelated records and global sort order are not certified by one lookup;
public immutable whole-file checksums retain that role.

After unchanged raw lookup caps and paired-span semantic merging, format-3 S
contexts are examined by chromosome, original position and original REF. These
already uniquely identify each merged context, so the fourth conceptual
`merged_ALT_list` tie never depends on physical offsets or omitted identifiers.
Singleton ALT token occurrence/order remains untouched, including repeated tokens
charged to implementation 7's work budget. An initial prototype incorrectly
canonicalized singleton ALT lists and changed the raw-token-work refusal; the
independent explicit fixture caught it. That code was removed before this frozen
candidate. Historical v1 artifacts are not final evidence.

The explicit format-3 receipt is
`reference_context_order="semantic-s-context-v1"`, plus `record_layout="typed-A7-S5"`
in reference manifests. Default format-2 behavior remains unchanged. The known
work-order qualification is classified, rather than called exact cap parity:

| Source order / cap | Frozen format 2 | Typed format 3 |
|---|---|---|
| Small quad first / 5 | ambiguous | ambiguous |
| Long noncompeting context first / 5 | context-work refusal | ambiguous |
| Either order / 200 | ambiguous | ambiguous |

All successful uncertainty outcomes have nullable canonical identity, effect,
EAF and orientation flags through R. No resource-exhausted proof publishes a
canonical target. The format-3 outcome is deterministic across both source orders.
This is a versioned resource-order difference, not a change to the scientific
meaning of ambiguity or context-required deferral.

## Verification actually completed

The API RED was observed against frozen implementation 7 before implementation:
five tests failed because `reference_format` did not exist; the old default-format
positive assertion passed. The final source and installed-library gates are:

| Gate | Actual result |
|---|---|
| Complete existing science suite, default 2 | 65 blocks / 1,170 assertions |
| Same existing suite using opt-in 3 | 65 blocks / 1,170 assertions |
| Existing suite plus new format-3 tests, default 2 | 72 blocks / 1,235 assertions |
| Installed registered public format-3 tests | 7 blocks / 65 assertions |
| Independent finite SNP oracle, all representations | 16,128 queries, zero violations |
| Independent whole-haplotype/contract oracle | 79,904 queries, zero violations; 2,104 resolved payloads independently valid |
| Independent explicit native contract controls | 297 rows: 156 resolved, 132 nullable, 9 fail-closed refusals |
| Mixed broad functional matrix | 190,392 complete native rows byte-exact across frozen7/reader2/typed3 |
| Writer source-order/sort-buffer matrix | 6 orders × 6 buffers: all six native format-2 payloads byte-exact; typed3 raw key multiplicities exact |
| Native structural controls | Existing core/block/projection plus new typed-format harness pass |
| Additional byte-corruption controls | 37 corrupt/compatibility cases pass, plus UTF-8/page/cap positives |

The opt-in-3 R science rerun uses copied original tests with the format default
injected to 3. Only **one metadata literal**, in
`test-source-reference-contract.R`, is changed from format 2 to 3. Scientific
expected statuses, coordinates, alleles, beta, EAF, raw source, counts and flags
are unchanged. The original source tests are unchanged. Both complete runs include
`test-equivalent-indel-context.R`; earlier small-run counts 64/1,163 excluded it
and are superseded by the complete final runs.

The entire C++ suffix from `Reference::harmonise` through scientific closure,
ALT–ALT/background guards, biallelic checks, trusted contract and projection is
**byte-identical** to scientific 7 (SHA-256
`47c1a1261aaff1421293d41089df74fedba37f264ed71ae78e930a9030b90d47`).
No canonical result, retained sign, EAF, ID or candidate count changes in the broad
matrix. The independent oracles check mathematics directly, rather than treating
old code as the scientific truth.

## Frozen identities and evidence

`results/mss-context-format3-stage7/source-pins.json` contains all source, native
driver, compiler and installed DLL hashes. Principal identities:

| Artifact | SHA-256 |
|---|---|
| `src/preservation_core.h` | `9c06526ec9185940d36e3002ff5d9c06a395d34e6a29d1cddb4b7ea119480011` |
| `R/preservation.R` | `60e1936f2e5378b534ed387351412cb46eff01b58ed6f18efd950f2f83a4dde1` |
| Installed prototype DLL | `d4f43e62b77a2135f94084ac3fe15105cbc6809fed388d2f4e7bc4ee9bc926da` |

Relevant final artifacts: `standalone-default2-final.log`,
`standalone-format3-final.log`, `installed-api-v2.log`, `native-suite-v2.log`,
`native-small-summary-v2.json`, `oracle-format3-v2/`,
`mixed-parity-summary-v2.json`, `writer-parity-v2/summary.json`,
`work-order-v2/summary.json`, `scientific-byte-parity.json`.
Native unresolved results retain existing false flags/position0 internally; R
converts all canonical fields/flags to NA. Raw inputs are retained separately.

## Prepared bounded measurements

`measure.py`, `bench.R` and the FAI-only driver are pinned by
`harness-pins.json`. Tiny guard success/censor selfchecks and a five-base public
build/query experiment exercise the executable harness. The first sandboxed `ps`
selfcheck refused closed; the successful test runs with process-metadata access
outside that sandbox and signals only its own newly launched process group.
No preparation processes or frozen libraries were signalled or changed.

The completed actual-prefix experiment reused existing immutable inputs:

- First 1,000,000 source sites, 31,408,384 bytes, SHA-256
  `3136a1b6bb0dd7e6bba90410bacf47ac2e10eaab0d07e190a38a1be8ed066372`.
- Full plain reference FASTA, 3,339,671,278 bytes, SHA-256
  `d88cc10556244b91b044294067ad166d9512bec6bbc104fb5c9f2d480af302fe`.

The same frozen installed DLL builds both formats in fresh paths. Timings separate
R/library setup, public prepare, native FAI+index, hash/copy/verification remainder,
and external wall; an independent FAI-only scan is reported separately, without
subtracting an uncontrolled-cache run as if it were an exact stage. Retained bytes
include every reference artifact and source copy; external shared FASTA is labelled
separately. Process-tree RSS and physical footprint and logical/allocated disk
peaks are sampled at 0.1s; macOS `time -l` additionally records exact single-command
maxRSS. Sampled disk peaks are explicitly samples, not complete filesystem-event
accounting.

Whole public generic queries use deterministic 256/2,048/10,000-ALT sites, with one
original row and all distinct concrete ALTs preserved. They include reference
verification, complete main-pool selection and bounded scientific proof, with all
19 public normalized fields and the exact original source compared across formats.
The proof cap is deliberately `max(10000,32*n)`; it is not a default-cap performance
claim. A source site with 10,000 A records has 10,000 raw C entries and 10,000 per
primary/alias I key, not one deduplicated candidate. Actual caps were 10,000 / 65,536 / 320,000 respectively. The original harness
flag `default10000_cap_tested=False` means no separate default-cap experiment;
for 256 ALTs the measured cap itself equals the default. A supplemental fast
native format-3 functional probe at 10,000 ALTs gives the same complete nullable
result at caps 10,000 and 320,000; cap 9,999 refuses the raw candidate count.
Format-2 default-cap timing was not repeated.

Every heavy command has an 8GiB process-tree RSS/physical guard, 160GiB free-space
floor, total 4GiB owned experiment disk cap and explicit deadlines. Exceeding a
bound is an operator-censored result, never a natural OOM or a successful benchmark.
FAI/build/query setup and uncontrolled OS cache are recorded. The shared host has
concurrent reference preparation; paired command order is fixed and no cold-OS-cache
or isolated-hardware claim is made. No new data acquisition is involved.

Independent format-3 review, native Windows CI and complete-reference capacity
admission remain required before default or production adoption. Root retains the
separate final-scientific-7 whole-file 1M/10M/40M throughput benchmark; these typed
storage/high-ALT experiments cannot substitute for it.


## Actual matched measurements and interpretation

All 15 heavy guarded commands completed successfully. None was operator-censored
or a natural OOM. The maximum sampled process-tree guard reading was
4,568,011,696 bytes; maximum aggregate owned experiment storage was
1,374,061,266 bytes; minimum free disk was 498,848,583,680 bytes. Source/DLL/harness
pins were checked before execution. Measurements are one fixed-order matched pair
per workload, with active background reference preparations and uncontrolled OS
page cache; small timing differences do not establish a general CPU speedup.

### Actual first-million-site reference prefix

The prefix contains **1,000,000 physical S occurrences and 1,242,071 physical A
occurrences**. A final streaming record audit verifies every emitted typed field,
record occurrence order, current ID, canonical event and original REF/full-ALT
pairing against format 2, in addition to exact source/FAI/rejected files and raw key
multiplicities. Thus the storage comparison includes actual multiallelic expansion,
not an assumed one-A-per-source-row panel.

| Metric | Format 2 | Format 3 |
|---|---:|---:|
| Retained index artifacts, including exact source | 232,818,833 B | 198,423,556 B |
| Allocated filesystem bytes | 235,593,728 B | 198,664,192 B |
| Sampled peak total build work directory | 289,446,572 B | 257,008,054 B |
| Native FAI + index | 11.487 s | 10.903 s |
| Public prepare, including hash/copy/verification | 26.589 s | 25.996 s |
| Hash/copy/verification + R remainder | 15.102 s | 15.093 s |
| Complete external command wall | 30.610 s | 29.894 s |
| Exact process maxRSS | 990,478,336 B | 1,007,173,632 B |
| Sampled tree physical footprint | 877,922,152 B | 888,948,512 B |

Retained index storage decreases **14.77%**. The external shared FASTA is
3,339,671,278 bytes and is excluded from both index totals. The independent
FAI-only full-FASTA scan took 2.009 s; that different uncontrolled-cache execution
is not subtracted as an exact native index stage. Peak work-directory bytes include
pending output and sorter files. They are a measured aggregate working footprint,
not an exact standalone sorter-scratch maximum; filesystem samples can miss a
short-lived maximum. This prefix shows similar build time and **no observed RAM
improvement**. It is incomplete and nonrepresentative, and does not certify full
reference capacity or whole GWAS-file throughput.

### Full public generic high-ALT queries

Each fixture is one eight-base-REF site with n distinct concrete eight-base ALTs,
current ID plus three aliases (including UTF-8), all retained in the source and
exactly one S plus n A typed records. The full public R path opens and verifies
the reference, performs main-pool C/I mapping and complete bounded generic proof,
and retains the original source strings. All 19 normalized fields and source
fields agree byte-for-byte across formats. Every result is
`source_allele_context_required`, with canonical identity/effect/EAF/flags nullable;
this is a complete conservative deferral path, not a resolved-hit throughput claim.

| Raw ALTs | Proof cap | Retained bytes 2 → 3 | Query processing 2 → 3 | Full command wall 2 → 3 |
|---|---:|---:|---:|---:|
| 256 | 10,000 | 632,668 → 38,755 | 0.010 → 0.004 s | 3.775 → 3.768 s |
| 2,048 | 65,536 | 38,104,972 → 314,216 | 0.650 → 0.010 s | 4.469 → 3.795 s |
| 10,000 | 320,000 | 901,794,681 → 1,579,087 | 18.734 → 0.045 s | 24.066 → 3.956 s |

The 10,000-ALT case additionally has reference open/verification 1.346 → 0.003 s,
exact process maxRSS 3,044,147,200 → 1,004,896,256 bytes, and sampled process-tree
physical footprint **4,568,011,696 → 883,918,672 bytes**. The physical-footprint guard
matters: resident RSS alone substantially understates the older representation's
memory in this run. Both remained under the 8GiB guard. Package/setup memory is
included and dominates the compact run; CPU times near milliseconds have coarse
resolution, and no precise small-case speedup ratio is inferred.

In this deliberately high-ALT shape, old A records repeat an O(n)-sized complete
ALT list O(n) times and main-pool context merging repeatedly processes it. Typed
A records omit that payload, while complete S context is preserved. These actual
whole public queries demonstrate removal of that storage/query burden without
changing source data or scientific results. The synthetic result does not imply
comparable gains for common biallelic studies or the 200,000-study corpus.

Final quantitative artifacts are `performance-summary-final.json`,
`prefix-1M-v1/summary.json`, `highALT-whole-query-v1/summary.json`,
`record-parity-final.json`, `default-cap-format3-functional-summary.json`, per-run
metrics/guard receipts and the evidence index. Benchmark files, source, DLL and
prior evidence remain immutable. Root has been notified that the heavy slot is
released for the separate implementation-7 10M/40M whole-file stages.
