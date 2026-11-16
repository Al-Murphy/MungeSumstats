# Final MSS7 whole-file benchmark harness: version-2 recheck

**Verdict: GREEN for the sealed version-2 harness and profile below.** The three version-1 gate/oracle findings are corrected and independently reproduced as refusals. The coordinator may proceed with the explicitly authorized matched 1M stage. This review did not execute any 1M/10M/40M benchmark or grant larger-stage authorization. Version 1 and its RED report remain unchanged.

- Harness: `results/mss-final7-whole-file-v2/`
- Seal SHA256: **`b2d6dd7d7f945f072ff2572d9eb0d3f737a8caea08ddd1e4c564aa779fbe4493`**
- Profile: `profiles/scientific7-format2-v2.json`
- Profile SHA256: **`6251cc132236435156b4dff4c3c934d7fe7505f263a04745def33d0f89983cd3`**

## Focused independent verification

The reviewer ran version 2's small self-check: **14 tests pass**, including every normalized-field corruption and the new approval/row-width/uniqueness adversaries. A separate reviewer script then repeated the original counterexamples using the actual current profile and actual scientific/installed-suite receipts:

| Probe | Version-2 result |
| --- | --- |
| Duplicate identical reference SNP rows | Refused |
| Undeclared extra normalized field | Refused |
| Explicit RED/NOT-GREEN scientific review | Refused |
| Stale installed-suite production revision | Refused |
| Stale profile production revision with unchanged core hash | Refused |
| Exact current GREEN scientific review and passing installed-suite receipt | Accepted by the pure approval function |

The last probe is a private test fixture; it did not invoke `run-stage`, create benchmark outputs, or replace the coordinator's execution approval. All static seal hashes were rechecked after testing. Evidence is in `review-mss-final7-whole-file-harness-stage7/`, particularly `selfcheck-v2.log`, `v2/focused-results.json`, `recheck-v2.py`, and `v2-comparison.json`.

The approval fix checks the current report's explicit anchored GREEN verdict prefix containing the exact production revision, the pinned report SHA/core, and the installed suite's matching source revision/DLL. The verifier requires the exact 19-key row shape and string values before value comparisons, rejecting both extra and missing-width rows. Reference validation now requires one source row per locus and one unique primary panel ID per row, in addition to the existing single-event, biallelic and nonpalindromic requirements.

## Source, workload and denominator confirmation

Version 1 and version 2 have byte-exact production source pins, installed library pins, reference inventory/files, stage input pins and resource limits. Only later checkout bookkeeping HEAD differs; both profiles declare and pin scientific production **`5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`**, core **`bbde5436366e686eeb35a48853117387b43b65618ca8308a87819e4382c86ceb`**, candidate DLL **`c1c0a8b29479c33d8b4cfa4646e46e39c4da436686f4ca2da12a402251003d50`**. Baseline remains frozen implementation 4 with DLL **`0a65467b2ef80c1c131bcc5f82a90fa0940b82e0eae0a2a33f71bf6591d981bd`**. The separately reviewed scientific and installed/native-CI evidence applies to that exact candidate source; this harness change does not alter MSS science.

The new preflight confirms **437,458 source rows, loci, IDs and eligible panel entries**, zero multiallelic loci and zero IDs with multiple events. The benchmark repeats this unique biallelic nonpalindromic chr22 SNP panel with exact synthetic beta/EAF. Every source row and all 19 normalized fields must be checked; both source copies must match original bytes. Setup, reference verification, processing, verifier and total operator time remain distinct. The 1M stage is the fresh implementation-4-versus-7 comparison; 10M/40M are separately authorized candidate-only scaling stages after a passing matched 1M pair. No 40M-unique-variant or legacy-40M comparison is supported.

## Carried guard assessment and scope

The reviewed runner, memory monitor and case guard are unchanged. R cases and the verifier retain finite deadlines, an 8 GiB RSS/physical-footprint/conservative process-tree target, 160 GiB filesystem free-space floor and owned-output growth caps. Only owned subprocess groups can be signaled. Loaded DLL/implementation checks, append-only output/refusal behavior, independent correctness verification and post-run pin verification remain in place. No active reference preparation, source, installed library or reference payload was changed by this reviewer.

The version-1 report's measurement qualifications still apply: polling can miss a transient, OS caches are uncontrolled, host work can overlap, and a single matched run is a bounded observation rather than a precision performance study. Optional hardening of missing external-time fields and explicit previous-stage profile/case linkage remains nonblocking within the current append-only named-profile workflow. There is no unresolved finding affecting this profile's scientific equivalence, row denominator or owned-resource execution plan.

**Proceed with the matched 1M stage using this exact profile and a separately pinned coordinator execution approval.** Obtain real successful correctness/resource receipts before claiming performance, and retain explicit authorization gates for larger stages.
