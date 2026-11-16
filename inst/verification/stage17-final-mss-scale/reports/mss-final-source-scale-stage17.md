# Final-source MSS whole-file scale verification

The final preservation implementation `03d11e92b417add0183dea4fb77716f5502d5735`
completed the 1M, 10M and 40M row workloads with every normalized row passing the
19-field oracle and both retained source copies byte-identical. The final 10M
and 40M runs were each executed once after independent operator review; neither
was substituted with an older implementation's timing. This is a measured memory
and correctness result, not a claim that legacy `format_sumstats` runs faster.

| Workload | Input rows | Processing seconds | External R wall seconds | Peak R RSS GiB | All 19 fields correct |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1M | 1,000,000 | 24.444 | 29.010 | 1.518 | 1,000,000 |
| 10M | 10,000,000 | 224.500 | 229.370 | 1.858 | 10,000,000 |
| 40M | 40,000,000 | 811.171 | 815.070 | 2.254 | 40,000,000 |

The 40M invocation exited zero, its fixed-size completion receipt is `accepted`,
no `failed-stage.json` exists, and root reverified all 81 small source, library and
control pins. Its processing time was 811.171 seconds; R including package setup
and reference open took 815.07 seconds. Exhaustive independent output verification
took another 258.781 seconds; the complete operator took 1,075.812 seconds before
its final commit. The completed profile occupies 12,228,014,479 logical bytes,
including its earlier stages and reserved bookkeeping. These are different time
and storage scopes and must not be added to a claim about R processing alone.

The exact production C++ header is
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`;
the loaded installed DLL is
`6f23605e06ac0b35fd40268b08190fd55e9a09e2fccded2f8cd7b629a3cd49eb`.
A later checkout contains documentation/CI changes only in relation to the pinned
production R/src/DESCRIPTION/NAMESPACE paths. Package version is 1.21.4, R is 4.5.3,
R/data.table/BLAS/OMP threads are one, and batches contain 100,000 rows. Profile
SHA256 is `158a8a8650a5364e28dd0eaaf61cc808e4a0e86c91d6a8f0f49249ca68d809d1`.

The workload repeats 437,458 unique real chromosome-22 biallelic non-palindromic
SNP panel records, with synthetic exact beta and allele frequency. Forty million
rows are not forty million unique variants or a complete real study. In the 40M
run 20,123,068 rows matched and 19,876,932 required allele swap; all corresponding
sign/frequency fields passed. Indels, multiallelic contexts, synonym ambiguity,
liftover and difficult real sources have separate reviewed functionality suites;
their correctness is not timed by this repeated-SNP workload.

Fresh R processes/application handles were used; OS page cache was uncontrolled
and never flushed. These are single runs on a shared host, not cold-disk medians
or a public latency SLO. The current final-source 10M run took 224.5 processing
seconds, compared with the historical legacy 10M run's 142.25 seconds: there is
no whole-pipeline legacy runtime speedup claim. Historical legacy peak RSS was
about 16.94 GiB, versus 1.858 GiB for this final-source 10M run; APIs/output work
are not identical, so treat this as workload-scoped memory evidence. The earlier
legacy 40M run was manually censored after excessive physical memory; it is not a
successful runtime denominator or a natural out-of-memory result.

Monitoring sampled synchronous known child/process groups every 0.25 seconds,
including macOS physical footprint. The 40M profile admitted an 8 GiB sampled
memory target, 160 GiB filesystem floor, 11,250-second whole deadline and finite
output cap. It did not establish a hard instantaneous memory/CPU cap or a general
orphan-aware production supervisor. No later benchmark, full reference preparation,
full count, index or public deployment is authorized by these receipts.

Root acceptance records:
- 10M: `01da94b5ca04159dc15260bf0e7d0acda9eed829b545aa6173b390b5c3a824a6`.
- 40M: `921b19a2bb961fb254d3ebafe25dea96a34d20a48bbdcbdf3e988e3f8f96818e`.

Both complete production reference preparations failed in the older merged-history
parser and have no successful seal. The separately reviewed empty-history fix has
not yet been shown to match the exact offending production record. Complete
reference count/index and the full 1,000-result corpus remain outstanding.
