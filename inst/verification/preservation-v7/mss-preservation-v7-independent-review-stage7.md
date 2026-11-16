# MSS preservation 7: independent scientific review

**Verdict: GREEN for the documented scientific contract of frozen `5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`.** No important scientific finding remains unresolved within its declared eligible-reference/context model. The earlier unsafe identity, sign, context-intersection and shared-background examples now withhold canonical measurements under explicit outcomes. Verified controls retain correct identities, effects and frequencies. This is approval to proceed to the coordinator's native CI gate; it is not a claim of Windows validation, deployment, adoption, or universal haplotype reconstruction.

All earlier RED reports remain unchanged and retain their original frozen-source scope. Format-3 proposals/prototypes are separate and are not adopted or reviewed as part of this verdict. The required runtime remains portable CPU C++17/R with native reference format 2.

## Frozen source and scope

Reviewed commit: **`5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`**, parent **`f2932c4abfeb9217cb2a620d38e1f602375085c3`**. A new `git archive` export was compiled independently; the author's checkout and installed libraries were not used to generate the independent oracle outputs.

| File | SHA256 |
| --- | --- |
| `src/preservation_core.h` | `bbde5436366e686eeb35a48853117387b43b65618ca8308a87819e4382c86ceb` |
| `R/preservation.R` | `d7bd275af108b49a9c588387c2098289b03fe602062032f92aa998cdf9ef2190` |
| `R/projection.R` | `51287b65e40a33ce9b616fb19ad78fdc0bc478cd636f8522f6e2fee0edc516a5` |
| `inst/PRESERVATION_CONTRACT.md` | `2d262653b54d64acc4fb711220fa0c16323329353068f560b06d2a1d57d3fa26` |

Scientific identities are normalizer **`mss-preservation-7`**, reference format **`mss-preservation-2`**, and projection **`mss-projection-2`** with normalization identity 7. The contract explicitly excludes arbitrary unlinked overlapping contexts and cannot certify completeness of the supplied reference or study metadata.

## Scientific conclusions

The main coordinate/input-ID mapping pool and its count retain their earlier meaning. A separate fixed-point rejection universe starts from both full initial pools, expands all current primary IDs and exact canonical-event coordinate bridges, and preserves original context locations. Known chromosome filters apply before expansion as documented. Unrelated overlapping C events do not bridge new IDs merely because they share a coordinate. ID cycles terminate. Proof-only REF–ALT hypotheses and paired ALT–ALT contrasts can reject a main selection, but cannot promote a replacement mapping or inflate the main candidate count.

For generic complex measured alleles, a missing source position now yields `source_coordinates_missing`, including the unique-ID/biallelic MNV sign counterexample. Its two known-position controls still produce the correct opposite signs and EAF transformations. Two distinct single bases retain ID inference: their REF–ALT normalization cannot shift, so closed canonical event positions exhaust their eligible REF–ALT anchors. The proof now quantifies over the completed universe rather than only the main pool.

S evidence remains paired by its own original chromosome, position and REF span, with concrete REF validated against FASTA. The revised policy distinguishes actual competing interpretations from insufficient metadata. `ambiguous`/`alt_alt_contrast` are supported by exact hypotheses/witnesses; `source_allele_context_required` is an explicitly conservative deferral when an alternate contrast cannot be excluded. Its necessary-condition tests do not pretend to reconstruct a full background/CIGAR alignment:

- Complex observed allele-length difference must match a difference between two distinct concrete ALT lengths for an ALT–ALT explanation to be possible.
- Single-base observations require equal-length original ALT sequences differing at exactly one base with the observed pair or its complement. This necessary condition can be checked without assuming that common alternate backgrounds preserve genomic anchor or REF-span length.
- Unknown ALT sequence/length in a multi-entry eligible context is handled conservatively. Duplicate ALT strings and entries equal to original REF are not a second allele.
- A context that jointly reduces REF and every concrete ALT to one pure SNP has a proved anchor; a different known source position excludes it. This prevents padding from creating false foreign-position uncertainty.

These boundaries resolve the known non-REF background plus shifted-repeat example and the single-base contrast on a common inserted background. Ordinary noncompeting triallelic SNPs, biallelic known-position indels/MNVs, and length-incompatible ALT pairs remain supported. Complete trusted source REF/ALT/strand declarations retain their separate behavior and can resolve measurements that the generic path must defer.

All source rows, identifiers and measurement tokens remain available. An uncertain normalized row is not missing source data, and a context-required outcome is not a claim that multiple mappings were actually proved. Consumers must retain that distinction when selecting training labels or reporting verified normalized coverage.

## Independent verification

Evidence root: **`/private/tmp/gwas-harvest-stage2/review-mss-preservation-v7-stage7/`**.

**Finite SNP oracle — 16,128 exact cases.** `oracle.py` independently enumerates REF/ALT sets and forward/reverse-complement measurement interpretations across four source-coordinate/ID modes, three frequency types, and minimal/left/right/doubly padded reference representations. All cases pass: **4,608 independently correct resolved payloads** and **11,520 correctly empty uncertainty outputs**. This includes the representation-dependent cases that blocked implementation 6.

**Contract-aware haplotype oracle — 79,904 cases across 120 contexts.** `haplotype-contract-oracle.py` compares complete alternate chromosome strings and separately enumerates tiny ALT-pair flank witnesses. It independently checks the necessary-pair exclusions/deferrals, all four source orientations and known/unknown chromosome/position evidence. It found **zero contract violations or false canonical calls**:

| Outcome/evidence | Count |
| --- | ---: |
| Uniquely validated resolved payloads | 2,104 |
| Coordinate-required deferrals, independently verified | 39,752 |
| Necessary-pair context deferrals, independently verified | 2,024 |
| ALT–ALT diagnostics, each independently witnessed | 8,192 |
| Ambiguous outputs, each independently witnessed | 544 |
| Allele mismatch | 15,208 |
| Coordinate/ID conflict | 12,080 |

All 17,736 independently found ALT–ALT explanations remain withheld, including rows whose stronger missing-coordinate policy takes precedence. The deferral counts are intended conservative behavior, not all proved multiple mappings. This deliberately varied synthetic matrix is not an estimate of real GWAS coverage or production deferral rates.

**Explicit counterexample, negative and resource controls — 308 native queries.** `contract-gold.py` has 297 predetermined queries; `closure-boundaries-summary.json` adds 11. Together they verify 160 exact resolved payloads, 134 empty uncertainty payloads and 14 fail-closed resource refusals. They cover:

- Distinct known-position indel identities and the opposite-sign unique-event MNV; missing-position complex rows, including an absent ID, require coordinates.
- Main-intersection context loss, partial aliases, coordinate-only derived-ID closure, exact-event linked IDs, proof-only REF–ALT alternatives and no promotion of proof-only mappings.
- All six source orders of the closed-ID counterexample, with full native fields identical; cyclic aliases terminate and multiple IDs for one event do not invent biological ambiguity.
- Non-REF backgrounds, repeat shifts, inserted backgrounds, unknown/`.`/symbolic ALT context, duplicate tokens, separate original REF spans and pure-SNP padding/foreign-anchor exclusion.
- Correct EAF/MAF/unknown-frequency handling, known/unknown chromosome eligibility, and skipping an unrelated known chromosome before an S lookup that would exceed its raw cap.
- Separate current-ID, event and location fanout limits; examined S-context limits; raw ALT-token and signature-work limits. Higher reviewed limits recover the expected verified/nullable outcome rather than changing allele direction arbitrarily.

**Public R and files — 117 independent assertions.** `public.R` verifies every deferred canonical coordinate/allele/ID/effect/EAF/direction flag is NA; exact resident source frames; signed known-position controls; EAF behavior; trusted source contracts; and normalizer/reference/projection receipts. Streamed editions retain original source-artifact MD5 and all rows at several parser batch sizes, with the new statuses and version 7 manifests. A deliberately insufficient work budget refuses the file and leaves no published destination; the original input and reference remain available.

**Compatibility and integrity.** Fresh standalone tests pass **65 blocks / 1,170 assertions**. Six source orders × six native writer payloads are byte-identical to parent implementation 5, including a valid 300-alias row and sort buffer seven that force multi-level external merging. The new reader accepts all six old format-2 indexes. Unsupported format-3 markers, malformed offsets and a queried malformed record refuse. Native parsing validates queried records rather than promising a whole-file structural scan; full public-R digest verification remains the separate immutable-reference integrity contract. An initial reviewer corruption control targeted an unqueried S record and incorrectly expected native global validation; the control was corrected to a mandatory queried A record, without a product/source change.

The coordinator's independently owned installed-package result was also read from `results/mss-preservation-v7/full-suite/summary.json`: **240 blocks / 2,013 assertions, zero failures/errors, two warnings and two skipped blocks**, DLL SHA256 `c1c0a8b29479c33d8b4cfa4646e46e39c4da436686f4ca2da12a402251003d50`, matching this source revision/core pin. That is the coordinator's actual installed-package gate, not a test rerun performed by this reviewer. Native Windows/other-OS CI remains pending at this report's verdict.

## Resource-order qualification and remaining limits

The 16× work budget counts additional raw lookup entries, examined ALT candidates, signature positions and normalization attempts. IDs, events, original locations and examined S contexts also have individual bounds derived from candidate_limit. Source-line length, allele length, repeat length and filesystem/process resources remain separate operator concerns; this is not a byte-level or wall-time guarantee. Existing format-2 full-ALT duplication remains a scalability limitation, not a scientific gain claimed by version 7.

**Tight work budgets can make uncertainty versus resource refusal depend on reference occurrence order.** `work-order-qualification/` uses the same two valid S contexts in two source orders. At candidate_limit 5, a small exact quad-allelic witness first yields `ambiguous`, while a long padded context first exhausts signature work before that witness. At limit 200, both orders yield `ambiguous`. Neither order publishes canonical measurements at either limit. This is scientifically fail-closed and does not contradict the contract's bounded examined-work policy, so it is not a blocker for this source. It is a real qualification: do not claim universal source-order invariance of resource outcomes, or assume a future format's changed physical-offset order preserves them. A semantic context-order/prefetch policy needs separate evaluation before format-3 equivalence/adoption.

The verified model is the explicitly declared current-ID/exact-event/coordinate closure over supplied reference contexts. It does not discover unlinked overlapping haplotypes, certify catalog completeness, harmonize every specialist entity, or replace missing source metadata with assumptions. Progressive metadata deferral and exact source retention are part of the accepted behavior.

No author source checkout, installed library, active GRCh37/GRCh38 preparation, user repository, cloud state or CI was modified. Only new private review fixtures/evidence and this report were written. The archived production pins were rechecked before the verdict. **Proceed to native CI on this exact source; retain the old RED scopes and keep any format/runtime optimization separately frozen and reviewed.**
