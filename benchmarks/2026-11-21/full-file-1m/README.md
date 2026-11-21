# Matched 1M rows through public `harmonise_file()`

**November 21, 2026 · macOS 27.2 arm64 · R 4.5.3**

Across three unprofiled pairs, the median complete public call fell from **25.950 s to 9.221 s (2.81×)**. Whole-worker median wall time fell **2.17×** and median peak RSS fell **29.50%**. All six editions agreed on stable artifact bytes and RDS fields. An independent oracle verified every one of the 19 normalized fields on all 1,000,000 rows in the first candidate edition. Exact normalized-byte parity extends that check to the other five editions. [Runner summary](raw/fullcall-1m-01/diagnostic-1m/summary.json) · [Independent SNP oracle](raw/fullcall-1m-01/snp-oracle-1M.json) · [Independent result review](raw/runtime/MATCHED-1M-RESULT-REVIEW.md)

| Metric | Experimental baseline | Candidate | Change |
| --- | ---: | ---: | ---: |
| Complete public call, median | 25.950 s | 9.221 s | 2.8142× |
| Whole worker wall time, median | 31.236721 s | 14.415484 s | 2.1669× |
| Package setup, median | 4.034 s | 4.112 s | |
| Verified reference open, median | 0.206 s | 0.229 s | |
| Process peak RSS, median | 1,800,880,128 B | 1,269,547,008 B | 29.50% lower |
| Process peak RSS, maximum | 1,846,886,400 B | 1,278,492,672 B | |

The three public-call ratios were 2.8142×, 2.6918× and 2.8098× in baseline/candidate, candidate/baseline and baseline/candidate order. Each call used one thread, a fresh R process and a fresh verified handle on a shared Apple Silicon desktop with 48 GiB RAM and uncontrolled OS file cache. RSS comes from complete Darwin `wait4` high-water measurements. The public-call interval includes input snapshotting and hashing, parsing, harmonization, serialization, checksum generation and publication. [Run conditions](raw/fullcall-1m-01/admission-1m.txt) · [Case evidence](raw/fullcall-1m-01/diagnostic-1m/matched-1M)

## Scope and correctness

The baseline was the experimental Git commit `a86e8c94d753630eff8c6a695ed2f240d1123d0d`. The candidate was `bfa192c83bb1fe15043c781d3b330f0f71bafb70`. Both were normal source installs with the same R, clang, SDK and dependency files. The result covers all intervening upgrades as a group. Published Bioconductor release comparisons, individual patch effects and the legacy R allele-flip component are outside this public-call measurement. [Baseline build](raw/fullcall-1m-01/baseline-build.json) · [Candidate build](raw/fullcall-1m-01/candidate-build.json) · [Common closure](raw/fullcall-1m-01/dependency-inventory.json)

The input repeated 437,458 real chr22 biallelic, non-palindromic SNP events with synthetic effects and EAF. The six outputs had 562,542 `match` and 437,458 `swap` rows each. The runner compared all eight stable non-RDS artifacts byte-for-byte and the stable contents of all six RDS manifests. The external oracle examined both source copies and all normalized rows in the selected edition. The separately recorded installed smoke covered hard and declared fixtures in reference formats 2 and 3. [Artifact digests](raw/fullcall-1m-01/diagnostic-1m/matched-1M/artifact-digests.json) · [Stable RDS check](raw/fullcall-1m-01/diagnostic-1m/matched-1M/manifest-verification.json) · [Smoke review](raw/runtime/INSTALLED-SMOKE-RESULT-REVIEW.md)

This is the single-size diagnostic. The full 23-pair matrix, 10M and 40M, warm handles, Windows, and indel/multiallelic timings remain open. The runner reports `target_status: not_applicable` for the diagnostic. [Benchmark guard](raw/fullcall-1m-01/matched-1m.jsonl) · [Oracle guard](raw/fullcall-1m-01/snp-oracle-1m.jsonl)

## Reproduce and locate evidence

1. Build both pinned commits into separate role libraries from clean `git archive` source tarballs with `R CMD INSTALL --use-vanilla`. Use the same dependency library for both and match the actual R, SDK, compiler, Rcpp inputs and flags in the build receipts. Capture fresh installed tree, library and runtime pins for the replication.
2. Supply the pinned Ensembl 104 GRCh37 chr22 FASTA and ordered panel. The baseline public `prepare_reference()` used build `GRCh37`, version `Ensembl104-chr22-bounded-prefix-unique-nonpal-SNPs`, reference format `mss-preservation-2` and a batch size of 100,000. Compare all seven stable member hashes, record the newly created manifest and verify the explicit FASTA binding in both installations. The pinned input can be recreated from the ordered panel with `scripts/benchmark-scale.py generate --panel "$PANEL" --input "$INPUT" --rows 1000000 --receipt "$GENERATION_RECEIPT"` and checked against the original input hash. [Input and reference pins](raw/assets/matrix-assets.tsv) · [Accepted reference](raw/fullcall-1m-01/reference-accepted.json)
3. Derive local configuration and SNP profile files from the recorded inputs below. Supply local physical paths, actual build/runtime/reference hashes and a fresh capacity and contention note. Keep format 2, explicit verified FASTA binding, EAF, beta, tab, the pinned batch settings, three pairs and no profiling. Use an admitted serial run with bounded wall time, disk and memory.

   ```sh
   BENCH="$REPO/benchmarks/2026-11-19/performance-upgrade"
   ORACLE="$REPO/inst/verification/fai-safety-stage8/whole-file/harness/verify.py"
   "$PYTHON" -E -s -B "$BENCH/run.py" run "$CONFIG" "$EVIDENCE_DIR" --suite diagnostic --profile none
   "$PYTHON" -I -B -S "$ORACLE" --profile "$PROFILE" --stage 1M \
     --edition "mss7=$EVIDENCE_DIR/matched-1M/pair-01-candidate/edition" \
     --receipt "$ORACLE_RECEIPT"
   ```

The `mss7` oracle selector requests the expected scientific behavior. The original config and profile show each input and output path. [Diagnostic config](raw/fullcall-1m-01/diagnostic-1m.json) · [Input cardinality](raw/fullcall-1m-01/matched-1M.keys.json) · [SNP profile](raw/fullcall-1m-01/snp-profile.json) · [Frozen harness pins](raw/assets/frozen-controls.tsv)

| Original location | Packaged evidence |
| --- | --- |
| `/private/tmp/mss-fullfile-20261121/fullcall-1m-01` | `raw/fullcall-1m-01` has the runner, oracle, restoration and smoke receipts. The six measured worker directories hold the settings, process accounting, R metrics and namespace inventories. |
| `/private/tmp/mss-fullfile-20261121/builds-01` | `raw/builds-01` has the normal-install guard records, logs and source archive pins. |
| `/private/tmp/mss-fullfile-20261121/runtime` | `raw/runtime` has reviews and the exact local metadata adapter sources for provenance. |
| `/private/tmp/mss-fullfile-20261121/native58-measured` | `raw/native58-measured` has separate supervision fixture and scoped process-closure evidence. |

Original receipts retain their absolute paths. The package contains metadata and hashes for the large input, reference, installed libraries and published editions. [Exact original-to-package mapping](SOURCE-MAP.tsv) · [Package checksums](SHA256SUMS)
