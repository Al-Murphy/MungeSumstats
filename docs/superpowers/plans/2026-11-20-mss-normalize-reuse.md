# MSS source-anchor normalization reuse implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development.

**Goal:** Avoid repeated identical source-anchor normalization without changing harmonization decisions.
**Architecture:** Four row-local orientation results, exact anchor/operand match and bounded64KiB ownedpayload; fallback computes unchanged original result. Proof and lookup work stays in its original order.
**Tech Stack:** Portable C++17, existing standalone/native fixtures, no R dependency for component run.
**Spec:** docs/superpowers/specs/2026-11-20-mss-normalize-reuse-design.md

## Global Constraints

- Base737ad247; product ownership src/preservation_core.h only, plus task tests/diagnostic/docs. Original source/DLL and benchmark fixtures immutable.
- Exact scientific results and error/cap/work ordering in both formats. Ordinary CPU/Windows, no public interface or referenceformat changes.
- <=4rowlocal results/64KiB additionalownedpayload, uncached fallback, no interrow cache or proof shortcut.
- Test/benchmark source draft first; no compiler/binary/R/native/AST/test/version/help/installation/cloud/process inspection until exact root-admitted command. No new data downloads or full reference traversal.

## Review Focus

Competing alternate/canonical anchors; absentcoordinates; budgetboundaryerrors; alias/multiallelic original context; cachepayload fallback/lifetime. Each is explicitly pinned by parity tests.

## Task 1: Native parity/work test and measured row-local patch

Files: src/preservation_core.h; new tests/native/test_normalize_reuse.cpp and task timing/instrumentation source. Consume exact controlheader and existing source fixtures; produce same Reference::harmonise results with fewer actual repeated sourceanchor normalizations.

- [ ] Draft finite real C++ parity/work harness, exact instrumentation/compile/run/resource proposal and unchangedproductpins. <=32unitrows percase,<=24paritycases, syntheticreferences<=64KiB perfile/<=1MiBtotal; optionaltiming loop explicitlyseparate, initially<=1024iterations/case. No product code until genuine baseline work regression observed.
- [ ] Root and independent review, then one admitted baseline build/run. Functional output parity may pass; actual work-count reduction must fail on old/old. No fabricated count or missinghelper failure.
- [ ] Implement only exact bounded row-local reuse, preserving fallback and every original proof debit.
- [ ] Separately reviewed covering native/parity/error tests and source/lifetime review; relevant installed/Windows integration remains separately gated.
- [ ] Compiler-identical plain control/treatment driver, same reference/data/order/caps, >=3serialpairs alternatingorder, no counters in timings. Record raw CPU/wall/peakRSS and correctness separately. Current mixed driver must be staged without its archived shadowingheader. Reject regression or no useful benefit.
- [ ] Commit exact covered files and task/final independent review. Full installed/matrix measurement follows runtime/input recovery; do not claim it from component tests.
