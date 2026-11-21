# MSS flip-index reuse implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development.

**Goal:** Reduce legacy allele-flip copying and repeated masks with exact behavior parity.
**Architecture:** Cache row positions within each original branch and table state. Preserve all existing scientific transformations and public contracts.
**Tech Stack:** R, real data.table, testthat/base-R fixtures.
**Spec:** docs/superpowers/specs/2026-11-20-mss-flip-indices-design.md

## Global Constraints

- Base737ad247; modify only R/check_allele_flip.R plus task tests/diagnostic/docs. User/frozen qualified checkouts unchanged.
- Exact row/order/key/class/value, warnings/errors/logs and sign/frequency parity; no new dependencies; ordinary Windows R.
- Test/harness draft first. No R/compiler/native/test/import/version/help/provider execution until separate admitted command; original source/control budgets remain charged.
- No fake data.table, Python translation, unverified runtime claim, installation, unrelated refactor, public push, cleanup or deleted evidence.

## Review Focus

Partial NA keep predicate; both alleles match with A2 override; integer positions after drop; effect coercion/negativezero; optional logs/approximatefrequency conditions. Each is an explicit parity case.

## Task 1: Exact parity, work-count regression, and measured patch

Files: R/check_allele_flip.R; new tests/standalone/test_flip_indices.R and task-owned diagnostics. Consume the unchanged original function from qualified source as oracle; produce identical function signature/results with fewer full-width count subsets.

- [ ] Cover up to 25 unit cases with at most 64 rows and 64 columns each. Include independent partial-flip cases for an integer `tmp` and a factor `flipped`, and compare all conditions and caller mutations to the oracle. Freeze the sources and runtime before measurement. Check the baseline full-width materialization count and verify the larger timing fixtures separately.
- [ ] After root+independent execution review, run the one baseline regression and read complete output. Functional parity may pass on baseline; operation-count assertion must fail for genuine full-width count materialization.
- [ ] Implement only row-index reuse, retaining original drop predicate and all branch/coercion order.
- [ ] After separate source review, run covering component tests and relevant legacy tests where dependencies exist. Report missing integration coverage exactly.
- [ ] Run serial matched original/patched component timing with same data/thread/runtime, both narrow/wide and sparse/dense cases, unprofiled time separately from allocation/work counts. Use >=3pairs and raw values; keep only a demonstrated benefit. This is component evidence, not full-pipeline speedup.
- [ ] Commit exact covered files, task review and fresh final review; root integrates only accepted bytes.
