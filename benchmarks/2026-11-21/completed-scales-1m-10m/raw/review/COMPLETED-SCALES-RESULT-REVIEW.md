# Independent results review of the retained 1M and 10M benchmarks

**Both scales pass the prespecified processing and memory target.** The candidate processed the 1M file 2.8243× faster and the 10M file 2.8415× faster by ratio of median processing times. Median peak RSS fell 27.19% and 38.57%, respectively. The two independent oracles verified all 19 normalized fields on every row. The original 19-case cold matrix failed at the disk floor.

## Measurements

I recalculated the table from the original twelve worker metrics.tsv and process.json files and checked the values against scoped-summary.json. Each scale had three serial pairs ordered baseline–candidate, candidate–baseline, baseline–candidate. All twelve processes had distinct PIDs, exited zero and reported complete positive Darwin wait4 ru_maxrss in bytes. Each cell with two values lists baseline followed by candidate.

| Measurement | 1M rows | 10M rows |
| --- | ---: | ---: |
| Median processing seconds | 22.498 / 7.966 | 211.391 / 74.394 |
| Ratio of median processing seconds | 2.8243× | 2.8415× |
| Processing ratios in pair order | 2.8458×, 2.8731×, 2.8243× | 2.9468×, 2.8935×, 2.5428× |
| Median whole-worker wall seconds | 27.390 / 12.394 | 215.467 / 78.768 |
| Ratio of median whole-worker wall seconds | 2.2099× | 2.7355× |
| Median whole-worker user CPU seconds | 22.498 / 11.325 | 175.210 / 70.393 |
| Median whole-worker system CPU seconds | 4.274 / 1.014 | 39.844 / 7.954 |
| Median peak RSS in bytes | 1,761,607,680 / 1,282,588,672 | 2,339,536,896 / 1,437,106,176 |
| Maximum peak RSS in bytes | 1,861,419,008 / 1,296,596,992 | 2,392,981,504 / 1,506,115,584 |
| Reduction in median / maximum peak RSS | 27.19% / 30.34% | 38.57% / 37.06% |

The performance rule requires a processing ratio of at least 2 and median peak RSS nonregression. The separate 8 GiB cap passed for all six candidate observations. The largest was 1,506,115,584 bytes (1.403 GiB). The whole-worker wall measurements include R startup, package setup and opening the reference.

## Correctness and provenance

- On the independent oracle's pair-01 candidate edition at each scale, source rows, normalized rows, rows with all 19 fields correct, expected and observed resolved rows, and the nonmissing effect and EAF denominators each equaled the full row count. At 1M, 562,542 matched and 437,458 swapped. At 10M, 5,187,962 matched and 4,812,038 swapped. Both deltas maps were empty, both unaccounted counts were zero and both source copies matched the input exactly.
- I checked each oracle’s recorded file sizes, SHA256 and MD5 before and after execution against the plan, original input inventory and selected edition's artifact inventory. All four file records agreed at each scale. The profile, verifier source, stage, physical edition identity and command selector agreed. The oracle JSON matched its captured stdout. The six editions at each scale shared one digest per stable file, including normalized.tsv. Each scale had six distinct raw manifest.rds digests. The manifest-only R process exited zero at each scale after enforcing absent on-disk paths and comparing every field except the validated creation timestamp. This connects the all-row check on one edition to the five other editions at that scale.
- The experimental baseline source was `a86e8c94d753630eff8c6a695ed2f240d1123d0d` and the candidate was `bfa192c83bb1fe15043c781d3b330f0f71bafb70`. Build receipts, config, role-tree hashes, role-specific DLL checksums, setting paths and loaded package namespace paths agreed. All twelve workers reported one thread, R 4.5.3 on aarch64-apple-darwin20, the same BLAS, LAPACK and locale, and identical version/path records for the 84 other namespaces. The original and current verification metadata agree on the scientific assets, role trees, harness and R toolchain. The verifier also records equality of its opening and closing inventories. The comparison applies to these source commits. A Bioconductor release comparison and any later HEAD require separate evidence.
- I checked the SHA256 values of the 97 original small controls, the 13 source pins and the selected-case evidence. The retained-edition verifier and the two oracle phase receipts each record completion, raw exit zero, confirmed closure and a reaped leader. The independent root observations record outer exit zero and absence of the corresponding process groups and sessions. The verifier's pending-oracle flags reflect its earlier execution. The separate oracle result and subsequent root observation close both checks.

| Evidence control | SHA256 |
| --- | --- |
| Reviewed verifier source | `fb61736c23dcde4bfe16fedfa450c72a67d8803273536022193d4d39004e2848` |
| verify-01/RESULT.json | `c5926660776d81c77eeab4d4fc9ab74776a197871bb89bd806b892efe2f59789` |
| verify-01/scoped-summary.json | `2fc45340f860a58acf1deba01ad8ba760235f6d665cfb5b1d03c3bfc92ad4650` |
| verify-01/oracle-window-plan.json | `4ebbaf3375fe28d2f6f5b6ca24676a5015882bc3e081a18d3298656659bfee23` |
| verify-observation.json | `e5b98e3cbd3f89c73d2622bd1467876ec7e04e531a979640c73890378f9bac40` |
| oracles-01/RESULT.json | `e2793d466915d52eca1a09d4bb6769151dd75b261f050189cf1bc2738dee3584` |
| oracles-observation.json | `a10a6080ab53e9bd4050856206aa952bfadd913856e0e33469f8e9cb443bd29a` |
| Original cold-failure-observation.json | `39a96694ac60f09bfad5e662c24ec2759fc88823149b2e94d2e3bd1628ecf995` |

## Interpretation and remaining scope

The measured workload repeats 437,458 distinct real chr22 biallelic non-palindromic SNP events with synthetic effects and EAF. It uses the format-2 reference, fresh R processes and fresh reference handles on a shared Darwin arm64 desktop. Asset pinning preread files and the OS page cache was uncontrolled. Each ratio describes three original paired measurements. The separate earlier 1M diagnostic contributes zero measurements to this review. Earlier hard and declaration gold fixtures address correctness in additional contexts. Indel and multiallelic throughput are unmeasured here.

The parent matrix retained 13 completed workers, including a single 40M baseline. The 40M candidate was censored. The original failure observation records supervisor exit 1, raw exit −15, confirmed closure and failure to persist a terminal JSONL record at the floor. This companion validates only the original twelve 1M/10M observations. The 40M comparison, full matrix, warm stage and C9 remain open, as do current legacy comparisons, profiles, lookup, full functional validation, Windows and any formal successor. I found no material gap in the two-scale result. Review used retained small evidence and source controls, without a new timing run, R invocation or large-payload hash pass.
