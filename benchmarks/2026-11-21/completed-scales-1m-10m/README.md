# MungeSumstats retained 1M and 10M benchmarks

The candidate processed the 1M and 10M datasets **2.8243× and 2.8415× faster** by ratio of median processing times. Median peak RSS fell **27.19% and 38.57%**, respectively. Both scales met the prespecified threshold of at least 2× processing speedup with median peak RSS nonregression. All six candidate workers also passed the separate 8 GiB peak RSS cap. [Independent result review](raw/review/COMPLETED-SCALES-RESULT-REVIEW.md)

| Measurement | 1M rows | 10M rows |
| --- | ---: | ---: |
| Median processing seconds, baseline / candidate | 22.498 / 7.966 | 211.391 / 74.394 |
| Processing speedup, ratio of medians | 2.8243× | 2.8415× |
| Processing speedup, pairs 1 through 3 | 2.8458×, 2.8731×, 2.8243× | 2.9468×, 2.8935×, 2.5428× |
| Median whole-worker wall seconds, baseline / candidate | 27.390 / 12.394 | 215.467 / 78.768 |
| Whole-worker wall speedup, ratio of medians | 2.2099× | 2.7355× |
| Median peak RSS in bytes, baseline / candidate | 1,761,607,680 / 1,282,588,672 | 2,339,536,896 / 1,437,106,176 |
| Maximum peak RSS in bytes, baseline / candidate | 1,861,419,008 / 1,296,596,992 | 2,392,981,504 / 1,506,115,584 |
| Median / maximum peak RSS reduction | 27.19% / 30.34% | 38.57% / 37.06% |

Whole-worker wall time includes R startup, package setup and reference open. Peak RSS comes from Darwin `wait4 ru_maxrss` in bytes. The largest candidate peak was 1.403 GiB. [Summary TSV](SUMMARY.tsv) · [Verifier summary](raw/verification/verify-01/scoped-summary.json)

## Workload and comparison

The experimental baseline was `a86e8c94d753630eff8c6a695ed2f240d1123d0d`. The candidate was `bfa192c83bb1fe15043c781d3b330f0f71bafb70`. Each scale contributed three serial paired runs in baseline–candidate, candidate–baseline, baseline–candidate order. Each worker used one CPU thread, a fresh R 4.5.3 process and a fresh handle on the format-2 reference. The workers ran on a shared Darwin arm64 desktop. Asset pinning preread files and the OS page cache was uncontrolled. The earlier standalone 1M diagnostic contributed zero runs to these estimates. [Original matrix configuration](raw/original/mss-fullmatrix-20261121/prep-01/matrix.json)

Both datasets repeat 437,458 distinct real chr22 biallelic non-palindromic SNP events with synthetic BETA and EAF. The measurements cover this workload and these exact experimental commits. Indel and multiallelic throughput, a Bioconductor release comparison and later commits require separate measurements.

## Correctness

The independent oracle checked all 19 normalized fields on every row of the original pair-01 candidate edition at each scale. Both source copies matched the input exactly. All rows resolved, and the nonmissing effect and EAF denominators equaled the row count at each scale. Both delta maps were empty with zero unaccounted deltas. [1M oracle](raw/verification/oracles-01/snp-oracle-1M.json) · [10M oracle](raw/verification/oracles-01/snp-oracle-10M.json)

| Scale | All 19 fields correct | Matched | Swapped |
| --- | ---: | ---: | ---: |
| 1M | 1,000,000 | 562,542 | 437,458 |
| 10M | 10,000,000 | 5,187,962 | 4,812,038 |

Before and after each oracle, all four bound files had identical recorded sizes, SHA-256 and MD5. Stable edition files, including `normalized.tsv`, were byte-identical across all six editions at each scale. All six raw `manifest.rds` digests per scale were distinct. The manifest-only R check required an absent on-disk path field and equality of every field except the validated creation timestamp. [1M edition digests](raw/verification/verify-01/matched-1M/artifact-digests.json) · [10M edition digests](raw/verification/verify-01/matched-10M/artifact-digests.json) · [Oracle closure](raw/verification/oracles-observation.json)

## Open work and evidence

The original 19-case cold matrix hit its filesystem floor and retained 13 completed workers, including one 40M baseline. The 40M candidate was censored. The supervisor exited 1, the workload received SIGTERM (raw exit −15), and closure was confirmed. The terminal JSONL record could not be persisted at the floor. The 40M comparison, full cold matrix, warm admission and C9 remain open. Current legacy comparisons, profiles, lookup diagnostics, full functional validation, Windows and any formal successor are also pending. [Original failure observation](raw/original/mss-fullmatrix-20261121/cold-failure-observation.json)

The original controls, later verifier and oracle receipts, and reviewed source are byte-exact in `raw/`. The earlier verifier's pending-oracle flags reflect when it ran. The later oracle receipts and root observations close the two scoped checks. [Source map](SOURCE-MAP.tsv) maps original paths to copies and records external identities. Full datasets, editions, reference payloads and installed software remain external. Run `shasum -a 256 -c SHA256SUMS` from this directory to verify packet integrity. The original commands and physical-path inputs are recorded in the receipts.
