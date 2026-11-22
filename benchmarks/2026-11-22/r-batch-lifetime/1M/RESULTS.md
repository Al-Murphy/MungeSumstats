# R file-batch-lifetime pilot: scientifically validated, modest improvement

The R-only lifetime prototype showed a small consistent improvement in this 1M
fixture: median paired operation-plus-closing-GC elapsed reduction **2.02%**,
whole-worker wall **1.44%**, whole-worker CPU **1.57%**, and peak RSS **0.45%**
(4,730,880 bytes median paired reduction). All three pairs improved on these
metrics. This is useful evidence for a small low-risk cleanup, not a major MSS
acceleration result. The product prototype remains uncommitted and unpublished.

The original producer phase still exited 1 after 729.399 seconds because I omitted
the frozen verifier's defaults merge. Its full 269 cases / 3,077 assertions passed
without skip/fail/error, and all six measured R processes had exited 0. A separate
verification-only completion ended 0 at 2026-11-22T01:32:02.412708+00:00 after 67.547 s.
It did not reinstall MSS, rerun any timing worker, write original outputs or
create original measurements.json. The original failure and both phases' receipts
remain intact. Comparison began only after this completion's outer guard passed.

Correctness acceptance covers 56 full stable-file hash checks (8 files ×7 editions),
stable RDS fields, all 19 independent SNP-oracle fields for 1,000,000 rows, original
full qualification, original 757 old/closed receipt identities and hashes, physical
installed trees, source 4 + FIX1 source 19 and clean unsignaled closure. Completion
returned all 11 allowlisted metadata files (289,946 bytes), source pins exact.

A is qualified installed 8b; B is an ordinary fresh installation of the authentic 8b
archive plus only the approved R overlay 32f. Three fresh-process pairs ran AB / BA / AB,
with identical 100,000-row / 8 MiB batch settings, input, reference, flags and one-thread
settings. All artifact/oracle/RDS verification occurred outside timed workers.

## Median values and paired changes

Seconds except RSS. A/B columns are separate medians; the percentage column is
the median of the three paired percentage reductions, 100×(A−B)/A. Positive means
improvement; negative means regression.

| Metric | A median | B median | Median paired reduction |
|---|---:|---:|---:|
| Operation elapsed | 13.849 | 13.526 | +2.37% |
| Operation + closing GC elapsed | 14.491 | 14.205 | +2.02% |
| Closing GC elapsed | 0.645 | 0.679 | -5.41% |
| Whole-worker wall | 21.406 | 21.097 | +1.44% |
| Whole-worker CPU | 19.981 | 19.668 | +1.57% |
| Whole-worker peak RSS (MiB) | 1008.871 | 1004.246 | +0.45% |
| GC CPU during operation | 2.201 | 1.534 | +30.30% |
| GC CPU during operation + closing | 2.845 | 2.212 | +22.25% |
| GC CPU through pre-receipt | 5.627 | 4.970 | +11.77% |

## Individual paired reductions

| Pair/order | Operation | Operation+closing | Whole wall | Whole CPU | RSS | Closing GC |
|---|---:|---:|---:|---:|---:|---:|
| 1/AB | +3.06% | +2.69% | +2.74% | +2.88% | +0.45% | -5.27% |
| 2/BA | +2.37% | +2.02% | +1.44% | +1.57% | +0.42% | -5.41% |
| 3/AB | +2.33% | +1.97% | +0.24% | +0.68% | +0.71% | -5.76% |

Closing GC regressed in every pair (median+0.035s,5.41% slower). Its cost is included
in the operation+closing and whole-worker results above. GC CPU during operation
plus closing fell 22.25%; GC CPU through the pre-receipt snapshot fell 11.77%.
Thus the observed net improvement survives inclusion of the closing collection.
Whole-worker wait4 includes startup, reference loading and small receipt writes.
The last GC snapshot precedes final small receipt serialization; its scope is
explicit rather than being called an exact whole-worker GC total.

All five public GC time components (user/system/elapsed/child-user/child-system)
for seven phase deltas are retained, including zero denominators and any negative
changes, in paired-summary.json, paired-details.tsv and median-details.tsv. GC
collection event counts were not instrumented and are unavailable. The public
[gc.time documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/gc.time.html)
defines these five times and notes possible timer-resolution underestimation.

This three-pair replicated chr22 SNP fixture does not establish whole-genome,
real-study, indel/multiallelic performance or a statistical confidence interval.
The broader functional suite and full scientific comparisons provide correctness
evidence, not a broad throughput estimate. No extrapolated total-corpus speedup
or new resource charge is claimed. Root handles the small R2 evidence archive.
