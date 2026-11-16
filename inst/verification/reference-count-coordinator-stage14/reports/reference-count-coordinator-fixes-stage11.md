# Sealed count coordinator — focused corrections, pending independent review

**Decision: no admission.** All four demonstrated independent-review findings are
corrected in a new no-hardlinks clone and branch. Behavioral source is frozen at
`0674d09930be48541d5e5306a0e72ebcb1164fd8`, branch
`fix/sealed-coordinator-review-stage11`, under
`/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-fixes-stage11`.
Fresh independent correction review remains required. Root alone may later admit
a separately pinned full-count profile. No full count, production index, full
profile/admission, provider/resource action or new performance qualification ran.

The preserved original is `547c7d46004c38f6a1f545748027af911ef28542`.
Its independent RED report SHA256 is
`034b2b1e80089cb934f7112c26b50d001424ae6856c6bd47a45e72fe8cad8313`;
review index SHA256 is
`7df58bbb6973abc31cb756b5aef9a1eada508cd0cd69c298f544c6e7e7e46ff2`.
Both review probe scripts were copied byte-for-byte into new evidence. Before any
fix, their logic reproduced all four findings on the unchanged original source
in this clone: publication-contained cases/profiles succeeded, oversized JSON was
content-opened before refusal, and the injected impossible native response wrote
success. The original proof directories were never reused. The explicit replay
harness relocates only the two hardcoded ROOT/CLONE string constants in memory via
AST; copied script bytes and all other logic remain unchanged. This limitation
and both literal replacements are logged, not represented as unmodified execution
against the reviewer's original output paths.
Byte-identical original coordinator/helper snapshots are pinned for auditing; the
historical replay receipts refer to paths in this new clone before source edits,
so compare their source pins with those snapshots or the unchanged original repo,
not with the later corrected files at the same paths.

## Corrections

1. The operator launcher now requires independently trusted expected GRCh37 and
   GRCh38 publication roots **outside parsed profile JSON**. Both CLI flags
   `--expected-grch37-publication` and `--expected-grch38-publication` are required;
   Python entry points require explicit `expected_publications`. There is no
   inferred-root fallback. Metadata-only placement checks reject profile, optional
   admission and case locations inside either expected publication **before any
   control content open or first mkdir**. Source aliases, symlink parents,
   traversal, overlapping roots and hardlinked external controls refuse. Parsed
   final roots must then equal the independent context before seal/control reads.
   The launcher binds/passes the same roots to the guarded private worker; internal
   worker/test entry points also fail closed without context. Placement refusal
   cannot mutate a publication even if a later guard/permission check would fail.

2. Expected and observed JSON sizes must be at most16MiB before content hashing
   or parsing. This applies to initial pinned JSON, parse reads, repeated producer
   control snapshots, and external profile/admission rehashes. Exact-cap JSON is
   valid. Reads are limited to the pre-open observed size, and opened-file/path
   identity plus full digest checks still reject growth, truncation or replacement.
   Two additional RED→GREEN regressions caught the repeated-hash variant after the
   initial size correction; those raw logs remain separate.

3. Native-response validation now applies necessary structural invariants before
   any planning receipt. Site/allele records have positive byte minima (typed S≥7,
   A≥13); key rows have a conservative minimum6 bytes. Zero populations require
   zero totals/maxima; positive maxima must leave minimum bytes for every other
   row. Source framing and maxima account for the26-byte header, five tabs per row
   and row separators. ALT disposition cannot exceed raw tokens/source bytes, and
   is exact when no positions are invalid. Empty ALT strings legitimately yield
   zero tokens, so no `raw_alt_tokens >= source_rows` assumption is imposed.
   Each record retains a key; each allele emits two raw coordinate keys. Raw/final
   byte differences and per-record fanout are bounded consistently with positive
   duplicate rows and maxima. Rejected bytes include the platform header and
   bounded row/status framing; nonempty FASTA requires positive FAI framing.
   The injected all-zero totals with positive records/maxima now refuse and write
   no worker success. This is malformed-derived-output testing, not evidence that
   the pinned real native algorithm miscounted. Scientific output is not filtered.

The original scientific/native code, original20 numeric fields, FAI bound,
I/P/both-order capacity formulas,8GiB/160GiB guard limits, finite deadline/output
caps,170GiB planning reserve and unconditional no-admission flags remain unchanged.
`full_reference_admitted=false`, `index_builds_executed=0` and
`decision=no_admission` still apply even when a planning envelope fits.

## Verification and scope

The frozen focused command used the existing environment, with process-metadata
permission for the unchanged guard:

```
python -m pytest test_coordinator.py test_coordinator_corrections.py \
  test_coordinator_integration.py -vv --basetemp=fix-evidence/frozen-focused
```

**105 passed:**44 existing coordinator cases,58 focused correction/boundary cases
and3 tiny guarded integrations. The initial correction suite observed35 failures
and12 passing positives before fixes, then47 passed. A focused compatibility run
exposed one missing-root exception-class mismatch (90 passed/1 failed), which now
refuses consistently. Repeated JSON rehash checks observed2 RED then2 GREEN.
Population/source/ALT neighbors observed4 RED then4 GREEN. Earlier intermediate
logs/fixture trees remain separate; no earlier output path was overwritten.

Placement tests trap every content open and compare source-tree names. They cover
both publications, profile/admission/case placement, source aliases, symlink parents,
hardlinked external controls, missing independent context and parsed-root mismatch.
Wrong roots may open only the top-level profile before refusing; source seals and
counted inputs remain unopened. Size tests distinguish expected/observed oversize,
booleans, exact-cap acceptance and post-validation growth before repeated hashing.
Malformed native neighbors cover zero bytes/counts, ALT contradictions, maxima,
remaining-row framing, rejection header/status bytes and FAI framing. Real native
positives retain empty/header-only input, all-rejected context, invalid positions,
empty ALT, indels, multiallelics, duplicate aliases, Unicode, minimal typed records,
CRLF and missing-final-newline boundaries.

Both newly owned actual tiny producer publications emitted all3 contexts each
(two original sites plus one historical-alias expansion). Count results match all
20 original fields from the exact original diagnostic on both builds. Two actual
small writer oracles independently match records/keys/offset/rejection bytes and
bound FAI bytes. Full publication filenames and file hashes remain unchanged by
the corrected count. Deadline/output censor tests retain failed guard evidence
without success summaries; only their owned tiny child trees are signaled.
No active preparation process is signaled. The final run creates exactly2 tiny
writer oracles; the earlier intermediate101-pass run created another2, recorded
separately. Count processes create no index output; full production builds remain0.

The unchanged C++ was not recompiled or re-sanitized in this pass. Known ordinary
counter/writer/original-counter binaries were copied with verified expected pins;
there was no new native concern requiring sanitizer replay. The earlier strict
sanitizer and FAI evidence retain their original source/scope. This is correctness
verification, not a new timing, memory-peak, Windows, prefix or whole-file
qualification. Root's temporary metadata timing hold was respected; later
correctness batches could overlap root's explicitly permitted PG/MSS correctness
work. Test elapsed times are not performance claims.

No full/active TSV, SQLite, FASTA, archive or staging content was opened, hashed,
scanned or changed. There was no broad task-base discovery in this pass. Access
used exact known source/control paths and newly owned fixtures only. No credentials,
network, provider, public/OAuth resource, environment/lock installation, production
MSS R/src/API/default, producer or guard changed. No .remember/tool memory was
read/copied, and HOME/home/CODEX_HOME remained inherited.

## Preservation, pins and evidence

The original repository remains clean at547c7d4. Original author/reviewer report,
index and source hashes were rechecked, as were all13 preliminary FAI artifact
pins. The original indexed fixtures and39-file diagnostic archive were not copied
or modified. The preliminary FAI scope correction remains intact: its historical
four logs contain76 tiny writer oracles, while count-process/full-production index
outputs remain0. Its old unqualified index field is not reinterpreted as zero test
writer invocations. Original author final ordinary/sanitizer totals42/40 remain
attached to their original runs, not this correction pass.

| Artifact | SHA256 |
| --- | --- |
| Corrected coordinator | 6d152cbd926f8949821a055f709023cecb8f20329b7f5873630b43c3b6fc069c |
| New correction tests | 84233531969fcb417f4ec3e25bd8d124c26d190d4d22b183813e2e5e57af6f54 |
| Unchanged native source | fb529bff15666905d39c5726d2f49638a15c8ee61d90880da1b13079a954a10b |
| Unchanged ordinary native binary | fa133e2c65bf9ed4ea0cb52791293b19b3a01fbead95f089ae13c0dff2b6c7fd |
| Original counter binary | fefe47feca68204d54915d302f5cd27daf21994b67d3dd4c965cc74ff778fe50 |
| Byte-preserved review probe | dc3281056bce56720725a71aa51c1eb85f582b3047270ab892c9cf51f458a6aa |
| Byte-preserved public probe | 8b5fa46138a08626e111471349d969b9c7aa8ad27222f4c90b9e2f94509160bc |

`coordinator-fixes-evidence-index.json` explicitly allowlists frozen source,
review-probe copies, observed RED/GREEN logs, selected baseline/final tiny receipts,
actual producer/writer/parity evidence and this report. It excludes itself,
symlinks, unrelated/active data, memory, and redundant unit fixture trees.
The owned16MiB+1 oversized proof and exact16MiB positive are losslessly gzip
compressed (16,328 and16,333 bytes); receipts pin original size/digest and verify
roundtrip decompression. Large sparse/repetitive proofs are not blindly copied
into the artifact allowlist. A final evidence-only commit may follow the frozen
behavioral source without changing its bytes. Independent correction review and
any later real full-count admission remain outstanding.
