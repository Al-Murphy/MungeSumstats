# Sealed reference count coordinator — independent review

**Decision: RED.** The coordinator must refuse publication-contained output cases
before creating anything. The reproduced placement and response-validation gaps
below also need resolution before considering a separately reviewed full-count
profile. No full count, production index, launcher or cost authority is granted.

Reviewed HEAD `547c7d46004c38f6a1f545748027af911ef28542`, behavioral source
`8989e74a782adc617853e0b8fdde7bc3b71b5797`, against original diagnostic
`1bf1a08685488798496e9d7978b1792dbd290a9d`. Review used a new no-hardlinks detached
clone at `/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11`.
It remains clean. The complete 475,759-byte original-to-HEAD diff was saved as
`complete-reviewed.diff` before either judgment. All source/document portions were
read; the complete index portion was parsed and all 1,188 explicit artifact paths
were read and verified, totaling 5,582,186 bytes. No broad task-base discovery or
blind author-fixture walk was performed.

## Findings

1. **P1 — The public launcher writes inside a sealed publication.**
   [sealed_count.py:333](/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11/sealed_count.py:333)
   validates `case.parent` but never compares the case with either final directory;
   `case.mkdir()` then mutates the publication. The internal tiny-fixture worker
   has the same omission at line281. In a new owned fixture, setting
   `case=<GRCh37 final>/count-case` completed through the unchanged real guard,
   wrote worker/count logs, guard and summary beneath the publication, and returned
   `counted_planning_only`. An initial sandbox refusal on `ps` had already created
   the directory inside the publication, so even an eventual refusal is too late.
   All original counted-file bytes stayed intact, but the publication tree was
   modified, contrary to read-only source handling and artifact placement rules.
   Evidence: `public-placement-probes-process-metadata.json`, the corresponding
   owned `public-case-in-publication-process-metadata/GRCh37/count-case/` receipts,
   and `independent-probes.json`. Check placement before the first mkdir, and apply
   the same rule to internal worker entry points used on tiny fixtures.

2. **P2 — External controls escape the publication placement check.**
   [sealed_count.py:317](/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11/sealed_count.py:317)
   loads the profile/admission into a separate `external` list. The publication
   containment check in `validate_seals` covers artifact, source-plan and producer-
   guard pins; it never checks these external paths after learning both roots.
   A correctly formed, pinned tiny profile at
   `<GRCh37 final>/operator-profile.json` was accepted by `load_profile` and a full
   public `launch`, which wrote a success summary. This is an unsupported exception
   to the promised controls-outside-publications rule. The current API learns roots by initially reading the trusted profile; this finding
   is the missing subsequent placement rejection, not a claim that reading a
   genuine top-level profile is itself counted-input access. A future launcher can
   accept independently trusted expected publication roots, check profile/admission/
   case placement before opening controls, then require parsed-root equality; no
   such mechanism exists in this candidate. The admission path
   has the same static omission; no synthetic full profile/admission was created
   to execute that branch. Evidence: `public-placement-probes-process-metadata.json`
   and owned `public-profile-in-publication-process-metadata/count-case/summary.json`.

3. **P2 — The control size bound is enforced after the entire hash read.**
   [sealed_count.py:111](/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11/sealed_count.py:111)
   calls `expected_pin` before checking the 16MiB JSON limit. `load_profile` uses
   the observed file size, so a huge file is fully streamed outside the worker
   guard before this limit can refuse it. A newly owned sparse 16,777,217-byte
   control produced a recorded `rb` content open and was then refused with
   `Control JSON exceeds finite bound`. This demonstrates ordering, not a large
   scan or an active-input access: only this 17MB owned fixture was used. Enforce
   the expected/observed size ceiling before content hashing, retaining stable
   identity/digest checks. Evidence: the `oversize_control_hashed_before_size_refusal`
   entries in both `public-placement-probes*.json`; no counted input was opened
   by this probe.

4. **P2 — Impossible numeric receipts pass native-response validation.**
   [sealed_count.py:260](/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11/sealed_count.py:260)
   checks types, rows, offsets and a few inequalities, but accepts four records
   with zero record bytes, three allele records with zero raw ALT tokens, and
   positive key/record maxima with zero key counts/bytes. It also accepts a zero-
   byte rejected file although the native contract always emits its header.
   The tiny native-output fault injection in `probe_review.py` preserves all
   fields and exact source/row identities, changes these totals to zero, and is
   accepted for both builds, writing `worker.json` and understated planning sizes.
   The real positive has R=64, U=170, K=89, N=10, J=11 and FAI upper=86; the corrupt
   response's corresponding values are all zero. Add necessary structural
   invariants before any planning receipt. This is a missing malformed-response
   defense using the same fault-injection technique as the author tests; it is
   **not** evidence that the frozen native algorithm emitted wrong counts, or a
   claim to sandbox malicious trusted Python. Evidence: `independent-probes.json`
   and owned `impossible-counts/case/` raw outputs and success worker receipt.

## Specification and quality assessment

The core implementation is otherwise well scoped. Both producer seals are
validated before the three counted inputs are content-opened. The normal gate
refuses a publication-contained artifact without opening the disguised input.
Exact reviewed source-plan/receipt semantics, producer/build/release/history,
emitted rows, and pre/post digest/stat identities are bound. Public counting,
alias accounting and counted-input pre/post hashes occur inside the unchanged
8GiB/160GiB process-tree guard. Both order formulas and the 170GiB planning reserve
match the brief; output/deadline censoring is failure. All successful probes still
report `full_reference_admitted=false` and zero count-process/production builds.

The native diff only adds the conservative FAI bound and another FASTA scan.
Name parsing matches the frozen private Fasta implementation; four maximum-width
uint64 fields plus native framing bound FAI serialization without exposing internals
or creating `.fai`. All 20 original numeric fields matched a newly compiled exact
baseline on both independently produced tiny publications. Production header
SHA256 remains `4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`.
No production R/src/API/default, original diagnostic source, guard or producer was
edited by this review.

The preliminary FAI evidence and all 13 pinned artifacts verified unchanged. The
additive correction appropriately distinguishes 76 preliminary tiny writer oracles
from zero count-process/full-production outputs. The author's final ordinary and
sanitizer runs separately account for 42 and 40 tiny writer oracles. The original
prefix timing remains attributed to its original source and excludes the new scan.

## Fresh verification and limits

- 16 targeted author cases passed on newly owned fixtures: positive read-only
  counting, seal/artifact early refusal, native failure/stderr/shape/row/integer
  checks, mutation during hashing/counting, overflow, exact boundaries and receipt cap.
- One targeted actual tiny producer/writer/guard integration passed, creating
  exactly 2 independent tiny writer oracles; each build emitted all 3 source contexts.
- 300 independently generated arithmetic pairs passed, including 1,200 both-order
  equality/one-byte-below reserve checks. Original 20-field parity passed on both
  actual tiny publications.
- Three public guard launches with process-metadata access completed: ordinary
  positive, misplaced profile, and misplaced output case. Initial `ps`-denied
  attempts are retained independently and were never overwritten.
- All author index pins verified. Author 88PASS logs and strict sanitizer evidence
  were checked without rerunning the broad suites: 81 successful calls have empty
  raw stderr and the single intended refusal has exactly its declared stderr.

These are tiny correctness and refusal probes, not a new memory/performance,
Windows or whole-file qualification. No active GRCh37/GRCh38 TSV/SQLite/FASTA,
archive or staging contents were opened, hashed, scanned or changed; no active
process was signaled. No new prefix, full count, full index, network, cloud,
provider, payment, credentials, environment installation or lock change occurred.
Both real preparations remain out of scope/unsealed for this review irrespective
of any later completion. Full-reference admission remains false.

Owned evidence is
`/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-independent-stage11/`.
Its `evidence-index.json` explicitly pins this report, review scripts, complete
diff, sources/binaries and owned tiny evidence. It excludes itself, symlinks and
unrelated data. Reproductions require fresh output paths and existing environment;
do not rerun scripts over their preserved evidence directories.
