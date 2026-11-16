# Sealed full-reference count coordinator — implementation, pending review

**Decision: no admission.** The operator-only coordinator is implemented and frozen
at `8989e74a782adc617853e0b8fdde7bc3b71b5797` on
`feat/sealed-reference-count-coordinator`, in
`/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-stage11`.
Fresh independent review remains required. No real full profile/admission was
created, no full-reference count was run, and no full production index was built.
Both production preparations were treated as unsealed throughout; their later
state, if changed, confers no authority in this task.

## Behavior and trust boundary

The local operator supplies an independently reviewed SHA256-pinned profile.
Expected source-plan hashes come from the reviewed acquisition/preparation context;
completed producer receipt/guard hashes come from that trusted completed context.
The coordinator does not calculate replacement expectations or authenticate JSON
received over HTTP. Full scope also requires a separately digest-pinned root
admission identifying the exact profile, frozen source commit, count-only action
and GREEN independent review. These local trusted records are not signatures.

Before any prepared TSV, derived FASTA or contig-alias contents are opened, both
builds must have real nonsymlink published directories, successful matching final
producer guards, exact reviewed source-plan/receipt fields, matching build/release,
preparation version4, complete same-release merged/withdrawn history, required
output identities and exact expected emitted rows. Pending directories, misplaced
source/artifact controls within a publication and malformed contexts refuse.
The complete manifest is compared, including reference mismatch policy. This is
important for the reviewed quarantine preparation; the counter does not change
normalizer predicates or replace unsupported/indel/multiallelic contexts.

The sole public CLI launches an owned worker under the unchanged process-tree
guard. Hashing every counted input, both native children, alias accounting and
post-hashing occur within that worker. The limits remain8GiB tree memory and
160GiB hard free floor, with explicit finite deadline/output cap. Counts are
sequential, with no row/variant/alias/ALT sampling, truncation or fanout reduction.
Censoring is failure and natural_oom remains false; it is not OOM evidence.
Worker/helper API calls used by unit tests are internal, tiny-fixture-only test
operations, not an unmonitored public CLI alternative.

Expected SHA256/size and stable device/inode/size/mtime/ctime/mode bind all three
counted inputs before and after counting. Coordinator, native source/binary,
production header, producer source, guard, memory harness, source plans and producer
receipts are also pinned before/after. Native failures, unexpected successful
stderr, truncated/invalid JSON, missing/extra fields, booleans/nonexact counts,
addressability/uint64 overflow, row mismatch and changed input/control pins refuse.
All native numeric fields and maxima are retained. Success summaries are exclusive,
appear only after guard success and final pin/resource checks, and are bounded to
the64KiB final-receipt allowance before writing. Failed cases remain separate.

## Logical planning

For each complete input, the unchanged formulas are:

```
I = T + R + K + 8*N + J + C
P = max(I, R + 2*U + J + C, R + U + K + 8*N + J + C)
A_then_B = max(P_A, I_A + P_B, I_A + I_B)
B_then_A = max(P_B, I_B + P_A, I_A + I_B)
```

C includes diagnostic `fai_upper_bytes`, exact quote-all/doubled-quote UTF8 alias
serialization with a conservative CRLF framing bound, the19-byte format marker,
and an explicit operator-supplied RDS manifest allowance/basis. The latter is a
planning assumption, not a measured or independently certified RDS size. Free
space is actually measured on the case filesystem after both existing inputs and
FASTA closures exist/authenticate. Allocation/directory/filesystem allowance,
concurrent reservations and170GiB planning reserve remain separate exact fields.
Both build orders and equality/one-byte-below reserve boundaries are tested.

Every apparent fit still reports `decision=no_admission`,
`full_reference_admitted=false`, `index_builds_executed=0`. These index fields refer
to count/production behavior, not the separate tiny writer test oracles. Any later
native build needs a separate independently reviewed build/memory profile,
maximum-line/sort-buffer assessment, fresh free measurement and post-first-build
remeasurement. No prefix-based small-file value or universal percentage is assumed.

## Fresh verification and observed limits

The final command was the existing environment's `python -m pytest -q
--basetemp=coordinator-evidence/final-pytest-v3`, with process-metadata access for the
unchanged guard. **88 tests passed**: the original29 native tests,12 FAI tests,
44 coordinator unit tests and3 actual guard/producer integration tests.
The final suite created42 tiny writer oracles:40 original/FAI comparisons and2
against real tiny producer publications. The actual producer fixture retains two
original source sites, the complete unsupported/multiallelic ALT lists and one
historical-alias expansion per build; all3 emitted rows are counted. Separate
native writer files match exact R/K/offset/rejection bytes, and actual FAI sizes
are bounded. No count process creates an index or `.fai`.

Observed TDD evidence is retained:32 missing-implementation failures then32 GREEN;
2 hardening failures (complete-plan binding and public unmonitored worker switch)
then39 GREEN;2 failures (source admission mismatch and invalid JSON refusal) then
41 GREEN;1 misplaced-artifact early-open failure then43 GREEN;1 final-receipt
allowance failure then44 GREEN. Existing-behavior checks include mutation during
hashing, native-time input mutation, receipt/source/build/version/history mismatch,
strict integers, randomized independent arithmetic, exclusive cases and source
immutability. Real guard deadline and output-cap tests censored worker plus native
child, retained raw receipts/logs and wrote no success summary.

The first integration attempt refused because sandbox policy denied `ps`; all
three resulting failures are preserved. Authorized process-metadata execution
then exposed an incorrect test expectation of4 emitted rows; actual producer
behavior correctly emitted3, so the test oracle was corrected. That failure is
also preserved. Later integration and complete suites passed. No sanitizer
failure was hidden or retried into a success claim.

Fresh C++17 ASan+UBSan builds replayed the original and FAI comparisons with
halt-on-error. All81 successful native calls had empty raw stderr; the one intended
format-suffix refusal had exactly its expected stderr. The replay creates40 tiny
writer oracles, preserves all82 stdout/stderr pairs, and includes arithmetic
overflow controls. Compiler flags were `-O1 -g -fsanitize=address,undefined
-fno-omit-frame-pointer` against the exact pinned production header. Native source
was not changed after the preliminary FAI slice; only new sanitizer binaries were
created. The prior source's2.864s prefix count is not a benchmark of this branch.

Tiny guard receipts contain incidental wall/memory samples. The final actual-tiny
success had only one sample (23,855,104 bytes conservative sampled tree memory);
its native child was too short-lived for a meaningful native peak qualification.
This is neither a maximum-memory proof nor a speed estimate for a prefix/full
reference. No new1M/full/timing qualification was run. Windows was not executed;
the unchanged operator guard supports macOS/Linux and fails closed elsewhere.
MSS Windows CPU portability and whole-file production memory/performance remain
separate qualifications. No environment, package lock, production R/src/API,
default format, original guard/producer, active data or original benchmark changed.

One initial overbroad `rg --files` discovery enumerated metadata-only filenames
under the task base, returning many historical fixture receipt paths and possibly
traversing active directories. It read no active TSV/SQLite/FASTA/archive contents,
performed no active file hashing, and signaled/modified no preparation process.
It confers no scientific or completion evidence. Subsequent access used exact
known source/control paths and newly owned tiny fixtures only. No further broad
active discovery or full-input action was performed.

## Preserved evidence correction and pins

`fai-slice-evidence.json` and all13 artifacts it pins remain byte-for-byte intact.
Its original unqualified `index_builds=0` needs the explicit scope correction in
`FAI-EVIDENCE-CORRECTION.md`: four preliminary logs executed76 tiny actual writer
oracles (12 initial RED +40 in the41-case GREEN +12 revised baseline RED +12 final
GREEN). That field means zero count-process index outputs/full production builds;
it does not mean no writer test ran. The final coordinator ordinary/sanitizer
runs separately add42 and40 tiny writer oracles. Earlier new-coordinator
exploratory reruns are preserved in their own logs and are outside these final
per-run totals. The original reviewed1bf1a086 diagnostic and39-file archive remain
unchanged and are not repackaged or reinterpreted.

| Artifact | SHA256 |
| --- | --- |
| Coordinator | ec226e236a5f71c8c297fc4076404795b95d25845d7a18a0f8225fa389b2e0c5 |
| Counter source | fb529bff15666905d39c5726d2f49638a15c8ee61d90880da1b13079a954a10b |
| Ordinary count binary | fa133e2c65bf9ed4ea0cb52791293b19b3a01fbead95f089ae13c0dff2b6c7fd |
| MSS header | 4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc |
| Producer module | c6327b05066a8a00508c724cd90ac44a0df2ccaa86da43831b43d888dd2b9fb5 |
| Unchanged guard | 1cf033bbbac586234310794fbd2bfd3a7f3bc4116a15b4296de870b052fe86be |
| Unchanged memory harness | aea36f87d4dd12ea41961f2cdd11198e527ccf17be2946a4351153cb96c0a0a5 |

`coordinator-evidence-index.json` is the explicit path/size/SHA256 allowlist for
tracked source, preserved preliminary FAI evidence, observed RED/GREEN logs,
final owned tiny fixture/publication/censor receipts and native sanitizer output.
Symlinks and unrelated task data are excluded. The report itself is pinned there;
the index excludes itself. No `.remember`, tool memory or active data is copied.
The final evidence-only commit may follow the frozen source commit without any
source-byte change. Root must assign fresh independent review before any further
real-input admission decision.
