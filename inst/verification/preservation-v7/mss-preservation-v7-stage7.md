# MSS preservation implementation 7: explicit source-context boundaries

Frozen source/test/recipe candidate: **`5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`**, parent **`f2932c4abfeb9217cb2a620d38e1f602375085c3`**, branch `fix/coordinate-contract-and-context-closure`, isolated clone `/private/tmp/gwas-harvest-stage2/mss-orientation-v7`.

This candidate addresses the scientific findings after implementation 6. It does not claim that a unique catalog event or matching short allele strings always establish an association contrast. Every successful call retains all source observations; insufficient context produces nullable normalized fields. Resource exhaustion still raises an explicit error and does not publish a partial file edition.

## Concrete failures and resulting behavior

**Missing-position complex alleles.** On FASTA CCAGG, one reference row at position 1, REF CCAG / ALT CCTG, normalizes to position 3 A:T. The same measured CCAG/CCTG strings yield beta +0.625 at known source position 1 and −0.625 at known position 2. A unique event ID cannot choose between these anchors/signs. Generic distinct complex DNA alleles now require a source position, even when the ID names one event. Missing-position rows return `source_coordinates_missing`, withholding canonical identity, effect, EAF and flags. Known-position controls and single-base SNP ID inference remain supported. Unsupported markers, same-allele observations and specialist preservation retain their established classification.

**Lost original contexts.** A source position 2 C/T observation can match a position 2 C:T reference event while a second position 1 AC / ALT GC,GT row with the same current ID supplies a genuine ALT–ALT C/T contrast on shared G background. Keeping only post-intersection canonical candidates loses that second original context. Ambiguity evidence now starts from both complete initial coordinate and input-ID pools, before their intersection.

**Linked proof-only REF–ALT events.** A source ID r may map directly to position 1 A:G. The same event's coordinate entry can reveal current ID q, whose ID mapping also includes position 3 A:G. For a position-missing A/G observation, the completed eligible universe contains two targets. Proof-only hypotheses now reject an otherwise selected mapping as ambiguous; they never supply a replacement mapping or change the main candidate count. Known source position 1 still resolves, and position 3 with source ID r still follows the original coordinate/ID conflict rule.

**Unmodelled shared backgrounds.** Exact joint trimming and individual reference-event normalization do not prove every composition of common alternate-background indels and repeat shifts. REF CAA / ALT GA,GAAA,CAAAA on CAAAA… can make source A/AAA at a shifted repeat position look like a reference insertion of two As, while the same strings compare deletion of one A versus insertion of one A on shared G background. An analogous SNP example is REF AA / ALT GCA,GCG,AG: source A/G can live on an inserted GC background. The new conservative `source_allele_context_required` status explicitly defers such necessary-only or unmodelled contexts. It does not call them proven ALT–ALT contrasts.

## Explicit scientific recipe

The shipped [preservation contract](../mss-orientation-v7/inst/PRESERVATION_CONTRACT.md) is part of the commit. It separates the mapping pool from the eligible rejection-proof universe:

1. Canonical coordinate/input-ID intersection, conflict behavior and `candidate_count` remain unchanged. The count means canonical candidates considered by these keys, not proof-graph size or catalog-wide absence.
2. Every current primary ID in the full initial pools participates in fixed-point ID closure. Exact canonical event coordinate bridges can find additional current IDs for that same event, including for ID-only input. Already loaded input-ID and initial coordinate lookups are reused. Known chromosome filters unrelated expansion.
3. Original S contexts remain separate by chromosome, original position and REF span. Source-order changes cannot pool unrelated ALT strings. Original REF must be concrete DNA and FASTA-valid. Only lightweight ID/event/context locators are retained across owned lookup batches.
4. Distinct closed-universe `(event, swap, complement)` interpretations or exact paired ALT–ALT witnesses can reject the selected mapping. Known-position hypotheses use the actual source anchor. Single-base position-missing hypotheses use canonical positions, which are complete because distinct single bases cannot trim or shift. Complex position-missing hypotheses are not guessed.
5. Necessary conditions distinguish uncertain context from proven ambiguity. For complex observations, two distinct concrete ALTs with the observed absolute length difference leave an ALT–ALT explanation possible. For single bases, equal-length ALT strings must differ at exactly one base with the observed pair or reverse complement; shared-flank signatures test this without inventing an anchor or positive REF-span length. Unknown ALT sequence/length also defers when a FASTA-valid context has at least two distinct non-REF ALT entries, including a genuine missing `.` entry alongside another ALT. Duplicate ALT strings do not create a second allele.
6. A fully concrete context that jointly reduces REF and **all** ALTs to one pure SNP has one proved genomic anchor. A known different source position excludes it, including equivalent common left/right padding. This exception does not strip a shared alternate inserted background. Ordinary noncompeting three-allelic SNPs remain resolvable.

Exact witnesses retain priority: `ambiguous` and `alt_alt_contrast` describe proved interpretations; `source_allele_context_required` describes conservative deferral; `source_coordinates_missing` describes the missing-position complex contract. All withhold canonical identity/effect/EAF/flags. Trusted complete REF/ALT/strand normalization is unchanged and can resolve an otherwise deferred generic observation.

This universe deliberately excludes overlapping contexts with unlinked IDs and no eligible coordinate/exact-event bridge. It is not an interval index or a claim about every possible global haplotype. Original provenance and trustworthy source REF/ALT/strand metadata remain preferred evidence. Indels, multiallelic sites and specialist observations remain in the source regardless of normalization confidence.

## Bounds and format compatibility

Raw lookup `candidate_limit` remains enforced. Closure admits at most that many unique current IDs, canonical events, context locations and examined S REF contexts. Additional work has a documented **16 × candidate_limit** bound, counting raw records before lookup deduplication, ALT entries, candidate signature positions and normalization attempts. Exhaustion raises a context-closure error; evidence is never silently skipped to publish a selected effect. IDs and events are deduplicated, including cycles; there is no ALT-pair Cartesian product or unbounded alias graph.

The work bound is not a byte/time bound for arbitrary allele strings or reference repeats. The existing line/batch/cache limits and operator memory/disk/time guards still matter. An external tiny functional probe verifies that a three-candidate/two-context case can exhaust the extra signature-work budget at cap 3 and fail closed, while cap 1,000 completes the exclusion proof and retains the same source and canonical match.

Normalization identity becomes **`mss-preservation-7`** in memory, file manifests/metadata and projection-normalization receipts. Reference format remains **`mss-preservation-2`**; projection algorithm remains **`mss-projection-2`**; parser identity and default writer remain unchanged. No reference reformatting, library mutation, cloud action or reference-preparation interruption was performed by this author.

| Frozen file | SHA256 |
|---|---|
| `src/preservation_core.h` | `bbde5436366e686eeb35a48853117387b43b65618ca8308a87819e4382c86ceb` |
| `R/preservation.R` | `d7bd275af108b49a9c588387c2098289b03fe602062032f92aa998cdf9ef2190` |
| `R/projection.R` | `51287b65e40a33ce9b616fb19ad78fdc0bc478cd636f8522f6e2fee0edc516a5` |
| `inst/PRESERVATION_CONTRACT.md` | `2d262653b54d64acc4fb711220fa0c16323329353068f560b06d2a1d57d3fa26` |

## Functional verification and exact denominators

New contract/closure tests fail against frozen implementation 6 for the intended missing-coordinate, lost-context and insufficient-proof behavior. The composed alternate-background/repeat case also has a recorded RED against the earlier working code. An early test expected a symbolic observation to be `unsupported_allele`; the existing wrapper correctly types it as `specialist_preserved`. The corrected control explicitly verifies that unchanged behavior.

Final portable standalone validation passes **65 blocks / 1,170 assertions**, including separately executed isolated testthat files. Native C++ gold also passes. Shared helpers live in the proper testthat helper file inherited from the parent test-only correction, avoiding the earlier cross-file-environment error.

Tests cover all four orientations, both known/unknown chromosome forms, positive candidate counts on deferred rows, unique-event opposite signs, exact known-position identity controls, partial historical aliases, coordinate-only current-ID recovery, ID-only canonical-coordinate bridges, all six source orders, rejection without proof-only promotion, chromosome exclusions, cycles, node/raw/work cap errors, byte-preserved streamed source, exact/potential background distinctions, trusted-contract positives and jointly padded pure-SNP negative controls.

Existing generic symbolic/third-allele cases now intentionally defer as context-required. Their trusted REF/ALT counterparts retain beta and existing conservative EAF behavior. Projection's trusted algorithm is unchanged; generic unresolved native rows do not acquire a projected effect. These changes are explicit scientific deferrals, not newly lost source observations.

The independent finite SNP oracle remains **16,128 complete public rows byte-exact** to mathematical gold across minimal, left-padded, right-padded and both-padded reference representations, including **4,608 resolved payloads** and **11,520 empty unresolved payloads**. All **24 native writer payloads** across the four representations are byte-identical to format-2 baseline.

A separate mixed native matrix compares every row/field against frozen implementation 5 (the earlier first implementation-6 candidate is exact on this matrix):

| Classification | Count |
|---|---:|
| Total functional queries | 190,392 |
| Complete result rows byte-exact | 180,890 |
| Retained resolved rows, all identity/effect/EAF/count/flags byte-exact | 4,293 |
| Prior resolutions → missing source coordinates | 110 |
| Prior resolutions → conservative context-required deferral | 120 |
| Prior resolutions → proved ambiguity | 24 |
| Prior unresolved `unmapped` / `allele_mismatch` → missing coordinates | 672 / 8,568 |
| New/changed deliberately small-cap refusals | 7 |
| Prior cap error → missing-coordinate unresolved row | 1 |
| New resolutions or retained-resolution sign/EAF/ID changes | **0** |

The seven refusal changes occur only at deliberately reduced limits 1–3, not the normal 10,000 limit; six formerly mismatched rows and one formerly swapped row now refuse incomplete proof. These are censored resource outcomes rather than row-level normalized successes. The classifier rejects any unaccounted change, candidate-count mutation on successful rows, newly resolved target or altered retained signed measurement. Native unresolved flags use their existing false defaults with position 0; R's nullable wrapper emits NA flags.

## Performance and adoption gates

No implementation-7 whole-file performance benchmark has run. The requested final-6 harness directory contains unfinished, explicitly gated preparation only; it has created no execution profile, authorization or run. Its code must be completed and copied to a new final-7 version directory with final reviewed source/DLL pins before use.

Earlier implementation-4/FASTA-offset and implementation-6 microbenchmarks keep their original identities. They must not be relabeled as implementation 7. Root is coordinating scientific review, a separate installed-package library/full suite, actual native Windows/macOS/Linux CI, and subsequent guarded matched whole-file benchmarking. No gate is inferred from the standalone pass.

Author evidence root: `results/mss-preservation-v7/`. Key files are `source-pins.json`, `standalone-summary-final-candidate.json`, `standalone-with-fallback-v5.log`, the RED/isolated GREEN logs, `native-gold-final.log`, `snp-padding-oracle/final-summary.json`, `mixed-matrix-summary-final.json` and full classified deltas, `writer-parity-final.json`, `work-budget-proof.log` and its reproducer, and `frozen4-5-integrity.json`. Original frozen implementation-4/5 source/DLL hashes remain unchanged.
