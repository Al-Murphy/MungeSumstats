# MSS performance implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox syntax for tracking.

**Goal:** materially improve measured end-to-end MSS time and memory while preserving scientific, parser and output behavior.

**Architecture:** independent legacy mapping, bounded native lookup and streaming/serialization changes in isolated worktrees, then serial integration and independent parity/performance review. Root coordinates shared interfaces and execution slots; the active reference build stays intact.

**Tech Stack:** R >=4.1, existing package dependencies, Rcpp, portable C++17, existing testthat/standalone harnesses.

**Spec:** docs/superpowers/specs/2026-11-19-mss-performance.md

## Global Constraints

- Preserve existing public APIs/schema/defaults and dependencies, ordinary CPU/R and native Windows support, R >=4.1 and C++17.
- Preserve synonyms/identifier fanout, indels, multiallelic/ALT–ALT ambiguity, sign/EAF orientation, conservative unresolved rows, source contracts and specialist preservation. Never fabricate declarations or weaken integrity.
- Preserve row order/duplicates, exact tokens, float64/NA/NaN/negative-zero semantics, quoting, Unicode/Windows text connections, batch membership, snapshots/checksums/receipts/audits and error/limit behavior.
- Native memory must be bounded independently of reference/study size. Do not change reference formats or require a rebuild or whole-index cache.
- Each agent edits only its owned worktree/files; root serializes integration and shared-interface changes. Do not modify the user's original source or untracked deep-profile artifacts.
- Initial source/control/Git admission is 128 MiB total under /private/tmp/mss-performance-20261119. No builds, package installation, R sessions, tests or benchmarks until a fresh root execution grant. Tests awaiting execution are explicitly pending; no passing baseline is assumed.
- Keep the original full GRCh37 build, limits and full Harvest/MSS objective intact. No providers, paid/public actions, pushes, main merges, cleanup or native retries.

## Review Focus

- Frequency-weighted contig selection and tied matrices must retain exact legacy behavior.
- Cached negative/rank lookups must preserve raw-entry limits, work consumption and corruption error order on warm calls.
- Scoped resident-page views must remain valid across eviction and boundary validation, including oversized UTF-8 keys.
- Typed parsing must retain decoded-byte batch grouping and audit order, including unselected columns and repeated EOF.
- Buffered encoding must preserve missing versus literal NA, negative zero, exponent thresholds, Unicode and platform text mode.

### Task 1: legacy chromosome mapping

**Worktree:** /private/tmp/mss-performance-20261119/legacy

**Files:** modify R/format_sumstats.R; create R/map_chr_style.R and tests/testthat/test-chromosome-mapping-performance.R. Do not edit preservation files.

**Interfaces:** internal map_chr_style(CHR, style) returns exactly the object from GenomeInfoDb::mapSeqlevels(as.character(CHR), style=style). The sole production consumer is the existing final chromosome assignment. No assembly argument or upstream filtering/liftover changes.

- [ ] Write oracle tests using expect_identical for NCBI/Ensembl/UCSC/dbSNP styles; duplicate names/order, factor and named input, empty/all-distinct input, autosomes/X/Y/MT, prefixed labels, unknown/scaffold/alternate/missing/empty labels.
- [ ] Add mixed-style frequency reversals and a tie. Assert equality with the full legacy expression and demonstrate that a naïve unique-only selection can differ.
- [ ] Request the root serial test slot; run the new test against baseline and retain the meaningful RED before implementation.
- [ ] Implement weighted complete-candidate ranking on distinct labels, restoring names from original chromosome values. Retain the exact legacy call for ties, zero candidates and compatibility cases from legacy-analysis.md. Use the helper only at the final assignment.
- [ ] Run focused parity and tiny available assembly/liftover integration cases under admission, retaining commands, results and limitations. Commit owned source/tests and self-review the diff.
- [ ] Report to the plan workspace task-1-report.md with status, commits, exact tests and pending evidence. Independent task review follows before acceptance.

### Task 2: bounded native reference lookup

**Worktree:** /private/tmp/mss-performance-20261119/lookup

**Files:** modify src/preservation_core.h; create role-specific tests/standalone/lookup_parity.cpp and tests/testthat/test-preservation-lookup-performance.R if needed. Do not edit R/preservation.R, src/preservation.cpp or generated bindings.

**Interfaces:** preserve Reference constructor, lookup(key,cap,work*) and harmonise(input,cap). An optional void Reference::prefetch_lookup_keys(const std::vector<std::string>& keys) may seed ranks only; communicate its exact memory/exception contract to root before the stream owner consumes it. Correctness may never depend on hints.

- [ ] Write baseline-versus-candidate complete-field/error/work oracle fixtures in both formats. Include first/last/missing keys, duplicate runs spanning page/median boundaries, changed caps and 0/N-1/N work across warm calls, negative lookups, corruption, oversized/shared-prefix UTF-8 keys, post-error reuse and page-cache-capacity-one eviction.
- [ ] Include aliases/current-ID bridges, repeat indels, symbolic and dense multi-ALT contexts, ALT–ALT evidence, positionless/source-declared queries and changed effects/frequencies/contracts under identical lookup keys. Do not cache Result or Candidate proofs.
- [ ] Request the root serial compile/test slot; retain a meaningful failing test or structural/counter regression against baseline before optimization.
- [ ] Implement portable scoped resident-byte comparisons, reusable page buffers and a lazy fixed pivot directory (maximum 8,191 slots of at most 64 bytes each; overflow goes to ordinary validated demand). Add an exact-key lower-bound rank cache bounded by both 1 MiB accounted storage and 4,096 entries including bookkeeping. Oversized keys are processed normally without retention. Preserve full format-3 row validation and raw occurrence/work/error ordering.
- [ ] If batch seeding is implemented, bound windows by 1,024 rows and 256 KiB retained hints with at most three hints per row; an oversized hint is declined before concatenation. Prefetch failures defer to ordinary demand and never become empty hits. No whole-batch Input/Candidate trees.
- [ ] Run focused dual-format/error/limit parity and instrumented tiny lookup controls under admission. Record counters separately from timing and preserve existing scientific kernels. Commit/self-review owned files.
- [ ] Report task-2-report.md with interfaces, byte accounting, tests/counters, commits and pending large-run evidence. Independent task review follows.

### Task 3: column-oriented parsing and exact buffered output

**Worktree:** /private/tmp/mss-performance-20261119/stream

**Files:** own R/preservation.R, src/preservation.cpp, any new local parser/serializer header, generated R/RcppExports.R and src/RcppExports.cpp, tests/standalone/registered.R and new tests/testthat/test-preservation-stream-performance.R. Also update only the new-header prerequisite in src/Makevars and src/Makevars.win so incremental Unix/Windows builds recompile preservation.o; no compiler-flag/default/toolchain changes. Do not edit src/preservation_core.h. Existing tests may be extended only where this track owns their behavior; keep new fixtures role-specific.

**Interfaces:** preserve public harmonise_file/harmonise_batch and current low-level parser APIs. New parser/serializer symbols are internal, with exact signatures in docs/superpowers/specs/2026-11-19-mss-stream-interfaces.md. Reference::harmonise remains sufficient; integrate optional prefetch only after root accepts the lookup track's concrete contract.

- [ ] Write old/new collector tests split at each byte and UTF-8/quote/CRLF/EOF boundaries. Cover identical row/byte grouping, empty/header-only/duplicate header, oversized row, wide unselected columns, absent mappings, ragged/NUL/quote/record-cap errors, truncated gzip/acquisition mismatch and no publication on failure.
- [ ] Write a baseline-writer byte oracle covering empty output, literal/missing NA, quotes/tabs/CR/LF/Unicode, negative zero, adjacent doubles, subnormal/max finite values, exponent thresholds, integer boundaries and logical missingness. Add permutations/duplicates, resident coercion/classes, formats 2/3, source contracts and specialist fixtures.
- [ ] Request the root serial Rcpp/test slot; retain meaningful RED before implementation.
- [ ] Share the syntax state machine with a bounded column sink, retain mapped columns while validating/counting every field, and keep at most one complete lookahead record. Preserve exact historical batch membership and raw snapshot writes. Reuse extracts and remove redundant single-batch binding/copies without changing resident coercion or source/projection checks.
- [ ] Implement bounded buffered encoding through R's existing text connection, with %.17g numeric semantics and exact legacy quoting/NA/text-mode behavior. Use a baseline-compatible fallback where platform parity is uncertain. Retain all integrity and per-batch mutation checks.
- [ ] Update registered symbols and standalone file inventories; run focused sourceCpp and registered-binding parity under admission. Commit/self-review owned changes and report task-3-report.md with exact tests, interfaces, commits and remaining evidence. Independent task review follows.

### Task 4: serial integration, full benchmarks and independent review

**Worktree:** /private/tmp/mss-performance-20261119/repo on feat/mss-speed-20261119.

**Files:** integration changes only through the owning implementer; root may add benchmark harness/docs under benchmarks/2026-11-19/performance-upgrade and verification reports. Root does not fix application code itself.

- [ ] Produce a task review package for each original task base a86e8c94..its final head. Have independent reviewers assess spec compliance and quality from complete diffs/evidence. Follow the bounded fix-loop workflow for real findings; do not repeat tests on unchanged hashes without a concrete concern.
- [ ] Integrate accepted commits serially. Reconcile optional prefetch API through the stream owner, then freeze source/DLL/compiler/runtime/input/reference/command pins for the combined candidate and exact baseline.
- [ ] Obtain fresh resource admission for package builds and complete parity suites. Run meaningful installed public API, standalone/registered, legacy, corruption, limits and platform checks. Retain Windows/native portability evidence for the new source; previous Windows results do not cover these changes.
- [ ] Obtain separate capacity/contended-run admission for full benchmarks. Run at least three fresh unprofiled 1M and 10M baseline/candidate pairs sequentially, plus one matched 40M pair and distinct/shuffled/coordinate/ID/negative/dense/mixed-variant workloads in both reference formats. Keep startup/reference-open, full processing, CPU user/system, peak RSS and allocation/counter diagnostics distinct. No OS-cache-flush claim.
- [ ] Compare full normalized/source artifact bytes, all normalized fields and stable receipt/audit/checksum fields plus independent hard-case scientific oracles. Record expected volatile path/time differences explicitly. Preserve every raw repetition and actual unique cardinality.
- [ ] Evaluate the >=2x median 1M/10M target with RSS nonregression and hard-case behavior. If missed, profile the measured remaining bottleneck and route a concrete additional change; report actual performance without relabeling scope. Full unique/genomewide coverage remains tied to complete reference availability.
- [ ] Dispatch one independent whole-branch scientific parity/performance review, using the most capable model. Address through its bounded fix wave, then deliver reviewed commits/reports into the user feature branch with scoped filesystem permission. No push/main merge/deployment. Do not mark the full Harvest goal complete on this subproject alone.
