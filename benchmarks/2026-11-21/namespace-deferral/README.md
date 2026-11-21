# MungeSumstats namespace deferral benchmark

In the matched 1M-row case, deferring the VariantAnnotation import reduced median package setup from **3.112 to 2.767 s (11.09%)**, whole-worker wall time from **11.055 to 10.582 s (4.28%)**, and OS peak RSS from **1,272,266,752 to 1,234,944,000 bytes (2.93%)**. Median processing changed from **7.117 to 7.018 s (1.39%)**, a 0.099 s difference in this run. [Independent result review](raw/mss-fullfile-20261121/runtime/NAMESPACE-BENCH-RESULT-REVIEW.md)

| Case | Rows | Median package setup, s | Median whole-worker wall, s | Median OS peak RSS, bytes | Median processing, s |
| --- | ---: | ---: | ---: | ---: | ---: |
| hard-cases-f2 | 32 | 3.091 → 2.743 | 3.707308 → 3.289779 | 1,016,905,728 → 943,079,424 | 0.004 → 0.005 |
| hard-cases-f3 | 32 | 3.073 → 2.756 | 3.686729 → 3.307992 | 1,031,061,504 → 961,118,208 | 0.005 → 0.005 |
| declared-f2 | 9 | 3.119 → 2.798 | 3.734650 → 3.346212 | 996,360,192 → 940,949,504 | 0.005 → 0.005 |
| declared-f3 | 9 | 3.115 → 2.795 | 3.736447 → 3.340149 | 1,044,840,448 → 933,920,768 | 0.005 → 0.004 |
| matched-1M | 1,000,000 | 3.112 → 2.767 | 11.054995 → 10.581757 | 1,272,266,752 → 1,234,944,000 | 7.117 → 7.018 |

Arrows show eager baseline → deferred candidate. The setup interval includes library-path selection, Rcpp and MungeSumstats loading, and the data.table thread limit. Peak RSS comes from Darwin `wait4 ru_maxrss`, in bytes. The tiny cases took about 4–5 ms to process and supply correctness and startup evidence. Their processing ratios are omitted. The diagnostic had no applicable formal targets. [Exact medians and units](SUMMARY.tsv) · [Worker source](raw/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R)

## Builds and measurement

The eager ANN-fixed baseline was `14ea7af51a19272c7792fbe582a3d42b291b4a1a`. The deferred candidate was `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b`. Both went through normal staged installation from clean Git archives into separate libraries under R 4.5.3, Apple clang 21.0.0 and MacOSX27.2.sdk with the same 96 common dependencies. All nine build phases closed and exited zero. [Build operator result](raw/mss-namespace-builds-20261121/operator-result.json) · [Build observation](raw/mss-namespace-builds-observation-20261121.json)

Each case used three serial pairs in baseline-first, candidate-first, baseline-first order. The roles had equal replication, with two baseline-first pairs and one candidate-first pair per case. All 30 fresh R processes used one thread and a verified reference handle, exited zero and recorded complete peak RSS. The shared Darwin arm64 host had no competing R or compiler processes observed at admission. Asset pinning preread inputs and the OS cache was uncontrolled. The exact 45 asset records and preparation, before and after installed inventories agreed. [Admission](raw/mss-namespace-benchmark-v2-20261121/admission.json) · [Benchmark observation](raw/mss-namespace-benchmark-v2-observation-20261121.json)

Each pair retained 85 baseline and 71 candidate namespaces. The same 14 names appeared only in the baseline, with zero candidate-only names. The roles shared 71 namespace versions and 70 paths. The one different shared path selected the role's MungeSumstats install. The first attempt stopped after two workers exited zero when the original 13-name gate saw the extra `otel`; its supervisor exited one. The installed DBI and testthat hooks explain why the earlier testthat reference probe included `otel` in both roles. The reviewed v2 wrapper adds exactly `otel` to the gate. The first attempt contributes zero observations to the medians. [V2 namespace policy](raw/mss-namespace-benchmark-v2-20261121/diagnostic/namespace-policy.json) · [First attempt](raw/mss-namespace-benchmark-failure-observation-20261121.json) · [Installed-hook audit](raw/mss-fullfile-20261121/runtime/OTEL-NAMESPACE-DELTA.md)

## Scientific checks

All six editions per case matched on stable output bytes. Stable RDS manifest comparisons passed on all five cases after checking the on-disk path field was absent and removing `created_utc`. The four tiny cases matched the literal oracles on all 19 fields in every edition. The independent oracle checked all 19 fields of all 1,000,000 rows on the first candidate edition, recording zero deltas, two exact source copies, 562,542 match and 437,458 swap statuses. The five bound input records agreed before and after the oracle. The 1M workload repeats real chr22 biallelic, non-palindromic SNP events across 437,458 distinct reference keys. The raw benchmark observation predates the separate successful oracle run. [1M oracle receipt](raw/mss-namespace-oracle-20261121/snp-oracle-1M.json) · [Oracle observation](raw/mss-namespace-oracle-observation-20261121.json)

The earlier unfiltered installed suite tested both source-aligned prototypes. Each passed 3,077 assertions in 269 cases across 89 active files, with zero failures, errors or skips, plus two matching existing warnings. Prototype checks covered first-use VCF loading, core semantics, parallel VCF and BP behavior. Source alignment binds those product files to the measured commits. The fresh builds provide tabular benchmark and parity evidence. [Functional review](raw/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-RUNTIME-REVIEW.md) · [Source alignment](raw/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-SOURCE-ALIGNMENT.json)

## Scope and packet integrity

The five-case diagnostic describes this exact commit pair and host. The earlier a86/BFA completed-scales performance cohort has separate estimates. C9, 40M, warm execution, fresh Windows qualification and publication remain open.

[SOURCE-MAP.tsv](SOURCE-MAP.tsv) gives original paths, packet paths, byte lengths and SHA-256 for all copied controls. External identities include the Git objects and archives, full 1M edition files, data and reference payloads, executables, DLLs and installed package trees. For those identities, receipt-recorded digests and lengths are carried forward when available. The tree hashes follow the retained inventory's tree-hash scheme. From the packet directory, run `shasum -a 256 -c SHA256SUMS` to verify every file covered by the checksum list.
