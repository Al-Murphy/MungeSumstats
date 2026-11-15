# Explicit reference contig aliases implementation plan

> **For agentic workers:** Use executing-plans for this small integration fix;
> a fresh reviewer checks it before integration.

**Goal:** Apply an authoritative assembly-report alias map to GWAS chromosome
labels while preserving original source fields and native reference coordinates.

**Architecture:** Add optional alias-to-FASTA named-vector configuration to
`prepare_reference()`. Store a checksummed portable table, reload it with the
reference, and map only declared input labels before the existing native call.
Production bioDB preparation emits canonical FASTA plus this explicit table.

**Tech Stack:** Existing R/C++17 preservation API; no new dependencies.

**Spec:** Accepted GWAS Harvest architecture and its C1/C7 production-reference
requirements; independent `production-reference-review-stage3.md` reproduction.
The user authorized implementation and MSS improvements. Label aliasing is an
implementation correction within that scope, not build inference or liftover.

## Constraints and review focus

- Original source fields, rows, measurements and bytes remain unchanged.
- Empty mapping preserves existing strict-name reference behavior.
- Many aliases may name one FASTA contig; each alias has exactly one destination.
- Unknown targets and remapping an existing FASTA contig to another are rejected.
- Alias table is hashed with reference files and carried into portable editions.
- Reference format2, normalizer4, projection2 and parser1 remain distinct from
  optional `mss-contigs-1` label processing; native code and native ABI are unchanged.
- Live-handle mutation detection and old reference manifests remain supported.

## Task1 — named input aliases and receipts

Files: `R/preservation.R`, `man/prepare_reference.Rd`,
`tests/testthat/test-contig-aliases.R`, `tests/standalone/preservation.R`.

- [x] Write failing fixtures for `1`, `chr1`, `NC_000001.10` mapping to canonical
  FASTA `1`, preserving raw CHR and all source/normalized measurements.
- [x] Test empty/unmapped labels, duplicate/unknown/shadowing alias rejection,
  many-to-one aliases, old references and live mapping mutation rejection.
- [x] Add `contig_aliases=character()` to reference construction; verify target
  names against the built FAI before atomic publication.
- [x] Persist/reload `contig_aliases.tsv`; map input chromosome strings only.
- [x] Emit checksummed `reference_contigs.tsv` and explicit mapping version in
  portable file editions and table provenance, preserving old inputs.
- [x] Run source/registered/native fixtures, installed suite and actual producer
  handoff; obtain independent review, then commit.
