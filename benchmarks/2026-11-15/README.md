These small receipts support review; large assets and raw logs remain in
`/private/tmp/gwas-harvest-stage2/results/mss`. Measurements use macOS/R4.5.3
and external `/usr/bin/time -l` RSS in bytes.

The matched panel contains 437,458 unique nonpalindromic biallelic real chr22
SNPs from a bounded Ensembl104 source prefix. Scale inputs repeat that panel
with explicitly synthetic BETA=.125/EAF=.25 and match/swap cycles. Adapted
legacy b4a2d21 uses exactly the same reference records; its QC is unchanged.
Added functionality is evaluated separately.

Successful unique-query runs measured 106.39s unoptimized versus 11.51s cached
end-to-end, and 102.297s versus 7.488s processing. The initial 1M diagnostic
that exited1 after output is excluded; formal runs use immutable snapshots.

| Rows | New cold process / RSS | Adapted legacy cold process / RSS |
| --- | --- | --- |
| 1M | 24.71s / 1,805,795,328bytes | 23.24s / 4,614,045,696bytes |
| 10M | 216.21s / 1,997,094,912bytes | 142.25s / 18,192,121,856bytes |
| 40M | 826.56s / 2,327,248,896bytes | Pending in raw matrix |

New has lower RSS but is slower than adapted legacy at 1M/10M. The frozen
cache-v1 40M result measured 2.17GiB RSS; strict identity revalidation is
tracked separately. It is not automatically a candidate-v3 scale claim.
Warm runs prime the complete identical input using the same verified handle.
External warm time/RSS includes priming; phase metrics separate it. OS caches
are uncontrolled. Tests, profiling and late verification scans overlapped the
40M run, so its timing is observational memory/functionality evidence.

Candidate-v3 preserves reference index version2 and identifies implementation
version3. It adds explicit trusted REF/ALT/strand and retains all equivalent
indel contexts before complementing EAF. Its installed suite passes 190 test
blocks/952 assertions, with zero failures/errors, two known warnings and two
missing-reference validation skips. Windows execution remains unmeasured.

For 14,428 actual PanUKB rows, trusted-forward candidate-v3 retains every row
and source byte, correctly normalizes all 14,422 FASTA-valid variants, and gets
all 12,686 eligible signed effects/EAF correct. Six reported REF mismatches
remain unresolved. There are 14,364 catalog calls and 58 sequence-only calls
without invented IDs. Historical aliases use dbSNP151 (2018), not a current map.

Synthetic legacy effect/EAF value-only correctness is 7/15 and 7/14; joint
canonical representation plus effect/EAF correctness is 5/15 and 5/14. New
achieves 15/15 effects, 14/14 EAF and 17/17 canonical calls and preserves 1e-400.
JSON receipts retain all denominators and separate unresolved observations.

Matrix receipts are settled-case snapshots. Historical legacy page-cache
labeling is corrected separately to zero without changing frozen raw receipts.
Strict verification binds new canonical key/ID/status to each source row and
checks the full legacy identity multiset. Corruption tests reject reordered
keys, wrong IDs and replacement by duplicated identities. Final pending outcomes
are maintained in the external task report.

Operational correction: the original RSS-only safety monitor was inadequate on
macOS. Legacy40M cold reached106,549,194,520bytes of compressed/private
`phys_footprint` while external peak RSS was22,648,717,312bytes. It was manually
censored, not a natural OOM or successful benchmark. Host swap was observed at
49,377MiB; no pre-run swap baseline was captured, so that entire amount cannot
be attributed to this case. Remaining legacy work is held until a tested guard
using max(RSS,physical footprint) is active. Small actual-allocation and complete
supervisor64MiB proofs are retained. Post-event elapsed timings are contaminated;
the new40M2.17GiB result remains explicitly a process-RSS observation, not a bound
on all host memory or filesystem cache use.
