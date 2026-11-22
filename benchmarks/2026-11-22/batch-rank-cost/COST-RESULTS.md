# Batch-rank cost diagnosis: preparation exceeds the saved lookup time

**Keep the selected16KiB implementation; do not integrate this prototype.** The
separate uninstrumented experiment was slower in all six native pairs, with median
native regressions7.414% matched and5.507% shuffled. This four-worker instrumented
diagnostic explains that outcome; it is not another speedup measurement.

| One1M-row worker per role/order | Matched | Shuffled |
| --- | ---: | ---: |
| Baseline ordinary-row stage |7.183196 s|10.664132 s|
| Candidate gathering/bookkeeping |0.089256 s|0.088727 s|
| Candidate prefetch |0.784772 s|0.780302 s|
| Candidate ordinary-row stage |6.722128 s|10.250274 s|
| Candidate clear/key destruction |0.002578 s|0.002603 s|
| Ordinary-row time saved |0.461068 s|0.413858 s|
| Total preparation/cleanup cost |0.876606 s|0.871632 s|
| Preparation minus row saving |+0.415538 s|+0.457774 s|

There is also a few-microsecond instrument bookkeeping residual in the total
clock. Timing is deliberately instrumented and one observation per role/order;
the six earlier uninstrumented pairs establish the performance decision.

## Search work improved, but that improvement was insufficient

Across each1M workload, total uncached midpoint probes fell22,214,082→10,906,614
(−50.90%), including3,501,867 prefetch probes. Ordinary-row probes alone fell to
7,404,747. Key-page requests fell29,160,321→17,820,144 (−38.89%), and offset-page
requests29,088,995→17,781,527 (−38.87%). These totals include preparation.

Key stream reads changed81,875→83,211 (+1.63%) in matched order and
743,335→696,013 (−6.37%) shuffled. Offset reads303→297 in both orders. Fewer
midpoint visits did not consistently mean fewer key-stream reads. Multiplying
these reads by16KiB keys/64KiB offsets describes application bytes requested,
including potentially short final pages; no physical/network-I/O claim follows.
Every candidate admitted all ten batch directories. Reported arena upper bound
was5,945,737 B; caller estimate9,600,000 B is not a strict heap bound. Matched
worker peak RSS was effectively flat; shuffled candidate RSS was about2MiB lower.

Gathering and clear stages had zero reference-counter events. All four stage
delta arrays reconcile exactly with their final counters. The candidate gains
2,000,000 exact rank hits, consistent with its initial C/I keys. It still performs
substantial fallback search and unchanged scientific resolution. The code's S
site-context lookup remains outside the directory; an S-heavy remainder is a
plausible inference for this biallelic panel, not a measured key-type breakdown.

The0.78s prefetch stage includes arena/copying, heapsort, deduplication and shared
traversal. The measurements DO NOT isolate sorting as the cause. Nor do they prove
that normalization, Candidate allocations or context closure now dominate the
remaining loop: directory binary searches and remaining C/I/S searches still
consume CPU. Algorithm changes must wait for actual hotspot attribution.

## Independent verification

Independently checked all17 returned files against exact size/SHA pins, exact
roster membership, source manifest and its four current files, all four full
13-field1M ordered gold digests, schedule/ordinal hashes, four distinct successful
wait4 receipts and commands, all counter sums and nonnegative stage clocks.
Reconciled each raw timing/process receipt with the nested outer proof; source
before/after outputs are byte-identical. Decimal arithmetic recomputes all values
in analysis.json. No scientific input body was read again locally.

The outer proof completed2026-11-22T07:28:53.412123 UTC with no signaled/unknown
cleanup, unchanged source/input identities and successful four-worker controls.
Proof SHA-256:
8c896c5ab8c0d6d1890e3d568cb88f5fd8a313d42ea7c778686666b76c514edc.
Diagnostic source manifest:
fd6638f3ddc4bed57a907f4313ed438387e4eba0175291420c512c440c4fc13b.
Baseline5881233f.../candidate9f71fa92... remain separate. No product integration,
new cloud work, compiler/native/R execution, or completion of the broader MSS
matrix is claimed by this analysis.

## Next focused profiling task — source recommendation only

Before implementing another optimization, use the existing native driver on a
future admitted cloud slot for two fresh5-second CPU samples: selected16KiB
baseline and experimental candidate, same matched1M panel and ordinary allocator.
Retain exact symbols/DLL or binary hashes, input/reference/gold and all13-field
digests. Use an available native sampling profiler; check its capability first,
rather than building a new guard framework. Keep samples separate from timing.

Attribute stacks within Reference::harmonise to: lower_bound_rank (ordinary cache
versus batch-directory search), lookup record decoding/merging, normalize/FASTA,
Candidate/set/map allocation/copying, and context_closure/alt_alt_context. Add only
test-only C/I/S lookup/probe counters if needed to identify the fallback remainder.
For the candidate's preparation, split one scalar per-call clock across arena/
key copy, heapsort, dedup and shared traversal; this isolates the measured0.78s
without per-row logging. Count preparation and ordinary work together.

This repeated-SNP control diagnoses the rejected experiment only. Before selecting
a product change, repeat the relevant attribution on both formats' existing
qualified mixed fixtures with explicit case counts, including shifted/declared
indels, multiallelic/ALT–ALT, aliases/multiple placements, missing coordinates and
negative queries. Retain every scientific path and cap/work/error-order control;
no SNP-only shortcut or skipped site context is proposed. Whole-genome/real-study
throughput and portable Windows performance remain separate acceptance work.

Actual compilation/sampling is deferred: original Pod compute cutoff applies,
and no further scientific work on it is authorized. All work here was small local
metadata reconciliation and source reading.
