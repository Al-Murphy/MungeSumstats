# Fixed-budget key-page experiment implementation plan

> For agentic workers: use superpowers:subagent-driven-development. Source, reviewed native counter execution, and promotion are separate gates.

**Goal:** Test whether smaller key pages reduce reference lookup work at the same bounded payload memory, before choosing another MSS processing optimization.
**Architecture:** Keep committed MSS f39ab1b unchanged. Build two standalone copies of its reference header against one identical parity/counter driver. Trial changes only the keys initializer from default64KiB×128 to16KiB×512; offsets, records, FASTA and all algorithms stay unchanged.
**Tech stack:** Existing portable C++17 header and lookup_parity.cpp, pinned Apple clang21/SDK27.2 for the local experiment; ordinary CPU/R and native Windows portability remain required for product promotion.
**Spec:** [Current profile diagnosis](CURRENT-MSS-PROFILE-DIAGNOSIS.md), SHA52b1c3f0062aae9fc7cd1dbd8be30571dd853feb4b6a26123857341f343d2336.

## Global constraints

- Preserve all source rows/tokens, float64 effects/frequencies and sign direction, indels/multiallelic contexts and one-to-many aliases; preserve ordered limits, errors and index formats.
- Keep key cache payload8MiB and total page payload28MiB perhandle. Additional container overhead must be disclosed and actual RSS measured before promotion.
- No GPU, external bcftools, WSL, new mandatory dependency or new monitor. Reuse existing qualified native guard for every compiler/scientific execution.
- All existing source, builds and accepted outputs remain unchanged. Scratch source/headers in `/private/tmp/mss-key-geometry-experiment-20261121`, no repository writes or publication for this experiment.
- Current profiles attribute514/3446 main-thread observations to lower-bound key-line stream reads; sampled ancestry proportions are not I/O-byte counts or predicted speedups.

## Review focus

- Long/UTF-8 keys crossing pages and the final short page: existing parity fixtures remain mandatory.
- Malformed typed records and precise error/limit ordering: existing format2/3 corruption and candidate/context fixtures remain mandatory; add new16KiB boundary characterization only if needed before promotion.
- Distinct sample and reproducible query sequence: select exactordinals and assert all32768 events unique, independent of repeated1M cycles.
- Float64 and source-token equality: compare every native Result field with double bits and the selected raw tokens, not just status or a lossy checksum.
- Performance interpretation: compare application page requests/read counts/requested bytes first; whole API timing and full scientific/platform qualification follow only if worthwhile.

## Task1: source-only matched driver and one-line trial patch

- [ ] Copy the existing standalone driver into scratch and add a bounded matched-reference mode described in KEY-GEOMETRY-BRIEF.md; preserve existing parity/pressure behavior.
- [ ] Write a trial patch of only the keys initializer. No production checkout mutation, no new public option.
- [ ] Source/syntax/control review, including exact output schema and actual input/variant selection. No contrived test asserting a configuration constant: scientific parity and measured counters test the behavior that matters.

## Task2: root-owned guarded counter experiment

- [ ] Capture exact baseline header from f39ab1b and stage trial with only reviewed patch; retain hashes and diff.
- [ ] Compile both roles with identical pinned toolchain/flags, assertions and MSS_LOOKUP_COUNTERS enabled.
- [ ] Run existing parity/pressure formats2/3 and require exact scientific parity/zero exits. Allocation-pressure and independent raw occurrence accounting remain real checks.
- [ ] Run both real distinct orders with fresh handles perrole/order, compare full outputs and selected ordering; report counters and requested key bytes at actual page widths.
- [ ] Independent result review; if unhelpful retain negative result and stop promotion. No claim of timing gain from counters.

## Task3: conditional promotion, still pending

- [ ] If mechanism supports it, commit the one-line candidate and meaningful test coverage in an isolated MSS branch; run matched unprofiled public API science/performance/RSS checks and current functional qualification.
- [ ] Require actual benefits and scientific parity before merge; native Windows speed claims need native Windows evidence. FullC9/40M/warm/publication remain separate.

Ruling: the user already authorized MSS acceleration, benchmarking, and portable C++ implementation. This small reversible configuration experiment is an implementation choice under that scope; no repeated consent or new deployment authority is needed. It adds a diagnostic mode rather than changing scientific behavior. Existing scientific characterization plus real counter measurements take precedence over a test that merely mirrors the cache constants.
