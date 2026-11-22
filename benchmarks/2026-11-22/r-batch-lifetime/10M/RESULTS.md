# 10M R-lifetime confirmation: modest CPU gain, slightly higher peak RSS

Median paired operation-plus-closing-GC elapsed decreased **3.63%**, whole-worker
wall decreased **3.47%**, and whole-worker CPU decreased **4.68%**. All three pairs
improved those timing metrics. Peak RSS **increased 0.40%**, with median paired
increase **4,308,992 bytes**; all three pairs had higher RSS. This supports a small
CPU/GC cleanup. It does not demonstrate a sustained memory improvement.

The trial ended exit 0 at 2026-11-22T02:30:32.924740+00:00 after 1071.431s, within its
fixed 1800s cap, with clean unsignaled closure. Three fresh AB / BA / AB pairs each
processed 10,000,000 rows of the repeated 437,458-site chr22 SNP fixture. No rebuild,
reinstall, resuite or repeated 1M benchmark occurred. A/B physical installed trees
and prior full 269-case / 3,077-assertion zero-skip/fail/error qualification were reused
and reauthenticated. Original 1M producer failure 729.399s and its later successful
verification remain preserved; this 10M phase is separately recorded as successful.

Scientific acceptance passed 56 stable-file comparisons (8 files ×7 editions),
stable RDS fields, one full independent 19-field 10M-row oracle with full-byte
transitive equality for all six new outputs, 789 protected old/closed receipt
identities/hashes, old 18 unretired payload identities, installed package trees,
namespace/options equality and source 4 + FIX1 source 19. There were no original output
writes or deletions. All 41 allowlisted metadata/log-tail files (403,415 bytes)
returned without absences; all 4 staged source pins and all 41 local file hashes/
sizes were verified. Root also independently verified these and archived metadata.

A remains installed tree 0ac1; B remains d732 with only the approved R-lifetime
source 32f overlay over the authentic 8b source. The measurement worker / public GC
instrumentation is unchanged. Both roles used identical 100,000-row / 8 MiB batches,
input/reference/flags and thread settings. One timed caller closing collection
follows operation return. Linux wait4 captures the entire fresh R process wall,
CPU and peak RSS. All full output/RDS/oracle inspection is outside these workers.

## Medians

Seconds except RSS. A/B values are separate medians; reductions are medians of
paired 100×(A−B)/A. Positive means improvement; negative means regression.

| Metric | A median | B median | Median paired reduction |
|---|---:|---:|---:|
| Operation elapsed | 135.456 | 130.444 | +3.63% |
| Operation + closing GC elapsed | 136.096 | 131.056 | +3.63% |
| Whole-worker wall | 143.074 | 137.922 | +3.47% |
| Whole-worker CPU | 131.483 | 125.347 | +4.68% |
| Whole-worker peak RSS (MiB) | 1023.176 | 1026.516 | -0.40% |
| Closing GC elapsed | 0.636 | 0.611 | +4.37% |
| GC CPU during operation | 20.439 | 14.390 | +29.61% |
| GC CPU during operation + closing | 21.077 | 14.996 | +28.84% |
| GC CPU through pre-receipt | 23.815 | 17.714 | +25.62% |

## Individual paired reductions

| Pair/order | Operation | Operation+closing | Whole wall | Whole CPU | RSS | Closing GC |
|---|---:|---:|---:|---:|---:|---:|
| 1/AB | +3.63% | +3.63% | +3.47% | +4.83% | -0.27% | +4.40% |
| 2/BA | +3.70% | +3.70% | +3.60% | +4.68% | -0.67% | +4.37% |
| 3/AB | +3.55% | +3.55% | +3.41% | +4.42% | -0.40% | +3.93% |

At 10M the closing collection was about 0.028s faster (median 4.375% reduction),
included in all totals. Operation-plus-closing GC CPU fell 28.84%; GC CPU through
the pre-receipt snapshot fell 25.62%. The median paired GC CPU decrease through
that snapshot (~6.101s) closely tracks the whole-worker CPU decrease (~6.148s).
These measurements are consistent with the intended lifetime cleanup; they do
not isolate the allocation producer. GC through pre-receipt excludes the final
small receipt serialization; whole-worker wait4 includes it.

All 47 metrics, including all five public GC time components for seven phase
deltas, all 141 signed comparisons and zero-baseline percentage denominators are
preserved in paired-summary.json, paired-details.tsv and median-details.tsv.
GC collection counts were not instrumented and remain unavailable. Public GC
elapsed timers can underestimate time because of resolution, as described by the
[official gc.time documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/gc.time.html).

The earlier 1M pilot showed 2.02% operation-plus-closing and 1.57% whole-CPU median
reductions with 0.45% lower RSS. The 10M result instead has a small RSS increase.
Keep these scale-specific results separate; do not multiply gains or claim a
major memory improvement. Both trials use the same replicated chr22 SNP panel,
not a whole-genome/real-study/indel/multiallelic performance sample. Three pairs
provide no formal confidence interval or corpus-wide throughput estimate. Broader
functional qualification supports correctness, not broad performance generality.
The product remains uncommitted and unpublished by this agent.

Actual guard growth was 14,112,683,520 workspace bytes and 303,104 root bytes within
the separate 20e9 / 1e9 limits. Peak total owned-tree RSS 1,539,538,944 bytes includes
external verifiers and is not the per-worker benchmark RSS reported above.
Original 186GiB/8GB/80GB floors and original 03:00/07:37/09:37 clocks remained intact.
All Pod workers and SSH returns from this agent are terminal; Root owns final
R2 report archival, GH issue updates and any product-integration decision.
