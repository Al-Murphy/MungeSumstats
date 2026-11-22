# Instrumented native batch-rank cost diagnosis — source only

No product semantic change, compiler/native/R execution or Pod call was performed.
Use only after Root's fresh fixed200s admission; END07:30 requires dispatch before
07:25:35 including5s kill+60s metadata return. Original compute/provider clocks
remain unchanged. This packet adds no guard/framework or resource authority.

Original benchmark was negative: all six native pairs slower, median+7.414%
matched/+5.507% shuffled. The selected16KiB product remains unchanged. These four
fresh instrumented workers diagnose cost, not a replacement speedup trial.

driver.cpp derives from the exact prior benchmark driver. Scientific row loop,
100k batches, seed/permutation/inverse, oneReference lifetime and13-field stdout
remain unchanged. Candidate gathering/limits/estimate/statistics match the prior
benchmark. Both roles now compile MSS_LOOKUP_COUNTERS; only candidate compiles
MSS_BATCH_RANK_BENCH. Baseline header remains5881233f...; candidate remains9f71fa92...
from the accepted prior source namespace. No header or input copied/modified.

Four per-batch clock stages are accumulated with scalar counters, no hot logging:
gather (includes caller estimate bookkeeping), prefetch, ordinary rows, and
directory clear plus key-vector/string destruction. Baseline records only its
ordinary stage; candidate includes all four. Total harmonise_ns includes snapshot/
bookkeeping gaps and is not required to equal the summed stage clocks. This
instrumentation changes cost; do not compare these times as uninstrumented gains.

Every stage records deltas in this exact order: rank_hits,pivot_hits,probe_entries,
key_requests,offset_requests,key_reads,offset_reads. Final cumulative values are
also emitted. Per-counter sum of four stage deltas must equal final values; all
counter activity occurs inside those stages. Gather/clear should have zero logical
reference events. Prefetch work counts too; never report just the row-loop saving.

These measure search/cache/ifstream operations. key_reads×16384 and
offset_reads×65536 are application bytes requested by this16KiB core, including a
potential short last page; not actual bytes returned or physical/network I/O.
No counters measure record/FASTA traffic, every string comparison, or system calls.
Heapsort/gather overhead is a hypothesis until the stage evidence is returned.

Future command under the existing accepted supervision:

    bash /workspace/mss-batch-rank-cost-diagnostic-20261122/run-native.sh

It creates fresh-only /workspace/mss-batch-rank-cost-results-20261122; compiles two
instrumented executables; runs matched baseline/candidate and shuffled baseline/
candidate once each. Existing meter/digest helpers retain wait4 and all1M rows.
All four digest JSONs must equal frozen panel.json (orderedSHA4abb49ca..., all13
field digests,562542match/437458swap). Expected order hashes remain matched
8137637191047640805 and shuffled9429961502938765523; each worker has10 batches and
1M rows. Candidate must report attempts/admitted10 and arena_peak<=32MiB.

Return16 small metadata files: source-before.txt,source-after.txt,compiler.txt,
binaries.sha256 and for each order-role its .timing.json,.process.json,.digest.json.
Leave the two binaries on Pod; scientific rows are streamed, not spooled.
No timing rerun, product integration, new scientific acceptance, cleanup or lifecycle
action is authorized by this source packet. Actual compile/output remains untested.
