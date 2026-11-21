# Matched 1M Darwin diagnostic result review

The first matched 1M installed public-call comparison passed. Median complete harmonise_file processing fell from **25.950 s to 9.221 s (2.814228×)**. Median whole worker time fell from **31.236721 s to 14.415484 s (2.166886×)**. Median process peak RSS fell **29.5041%**. The six editions passed byte and stable RDS parity, and the separate oracle verified all 19 output fields on every row of the selected edition.

## Measured comparison

The roles were fresh normal installs of baseline a86e8c94d753630eff8c6a695ed2f240d1123d0d and candidate bfa192c83bb1fe15043c781d3b330f0f71bafb70. All six workers had distinct PIDs, ran serially, exited zero, and reported one million source and normalized rows. The run used three unprofiled pairs on Darwin arm64 with format 2, verified FASTA binding and EAF/beta/tab settings.

| Pair | Observed order | Baseline processing | Candidate processing | Baseline ÷ candidate |
| --- | --- | ---: | ---: | ---: |
| 1 | Baseline, candidate | 25.950 s | 9.221 s | 2.814228× |
| 2 | Candidate, baseline | 26.695 s | 9.917 s | 2.691842× |
| 3 | Baseline, candidate | 21.796 s | 7.757 s | 2.809849× |

| Per-role statistic | Baseline | Candidate |
| --- | ---: | ---: |
| Median complete public-call processing | 25.950 s | 9.221 s |
| Median package setup | 4.034 s | 4.112 s |
| Median verified reference open | 0.206 s | 0.229 s |
| Median whole worker wall time | 31.236721 s | 14.415484 s |
| Median peak RSS | 1,800,880,128 B | 1,269,547,008 B |
| Maximum peak RSS | 1,846,886,400 B | 1,278,492,672 B |

The RSS medians and maxima use complete Darwin wait4 ru_maxrss records from all six workers. Processing, package setup, verified reference open and whole worker wall are the distinct recorded boundaries. The speedup is the ratio of the two processing medians. It describes the full a86-to-bfa revision range, including prior upgrades. [Runner summary](../fullcall-1m-01/diagnostic-1m/summary.json)

## Integrity and scientific verification

All workers reported R 4.5.3 on aarch64-apple-darwin20 with identical R engine, libR, BLAS, LAPACK and locale observations and one thread. Each worker reported its selected MSS DLL path and MD5. The review independently matched both actual DLL SHA-256 values to the configuration. The workers loaded the same 85 namespaces after normalizing the selected MSS role path, with 69 from the restored common library, 15 from the R framework and MSS from the selected build. The 20 registered asset pins and both complete 695-file MSS inventories agreed before and after the run. [Before pins](../fullcall-1m-01/diagnostic-1m/pins.json) · [After pins](../fullcall-1m-01/diagnostic-1m/post-pins.json)

Each of the six editions contained the same nine regular artifacts. All eight non-RDS artifacts matched byte-for-byte in the runner's comparison, with matching sizes, SHA-256 and MD5 in the per-edition inventory. Both source copies per edition matched the 46,106,766-byte pinned input. The six RDS files had distinct raw hashes. The successful RDS check compared their stable contents. The review checked each actual edition's inventory, artifact sizes, checksum table and row metadata against the receipts. [Artifact digests](../fullcall-1m-01/diagnostic-1m/matched-1M/artifact-digests.json) · [Stable RDS check](../fullcall-1m-01/diagnostic-1m/matched-1M/manifest-verification.json)

The independent oracle checked the first candidate edition against the repeated chr22 panel of 437,458 eligible biallelic non-palindromic SNP events. All 1,000,000 source and normalized rows and all 19 normalized fields passed, with 562,542 matches, 437,458 swaps and zero unaccounted deltas. Exact normalized byte parity extends the result to all six editions. The runner's embedded oracle count was zero because this independent check ran as its own phase. [Harness verification](../fullcall-1m-01/diagnostic-1m/matched-1M/verification.json) · [Independent oracle](../fullcall-1m-01/snp-oracle-1M.json)

## Bounds and scope

The benchmark guard completed with raw exit zero, confirmed group closure and a reaped leader. It recorded 136.098902 s elapsed, 133 group samples with zero misses, peak sampled group RSS of 1,878,081,536 B, 1,525,325,824 B peak allocated and 1,493,361,730 B peak logical disk, and a 290,301,534,208 B free-space minimum. The admitted limits were 1,800 s, 8 GiB group RSS, 3 GiB disk and a 186 GiB continuing floor. The separate oracle likewise completed with raw exit zero, confirmed closure and a reaped leader in 11.791142 s. It peaked at 1,575,534,592 B sampled group RSS within its 600 s and 4 GiB limits. [Benchmark guard](../fullcall-1m-01/matched-1m.jsonl) · [Oracle guard](../fullcall-1m-01/snp-oracle-1m.jsonl)

The benchmark used fresh processes and verified reference handles with an uncontrolled OS page cache on a shared desktop. The observed case exceeds 2× median processing improvement with lower median RSS. The diagnostic summary records target_status = not_applicable. The full 23-pair C9 matrix, 10M and 40M runs, separately scoped warm measurement, Windows, and indel/multiallelic timings remain outside this result. The comparison does not isolate a speedup for any individual patch. The review recomputed the metrics from all six persisted process receipts and cross-checked the independent oracle receipt and stdout. It ran no benchmark, R process, oracle or network request.
