# Keep 16 KiB pages; retain the independently qualified R cleanup

The extra 8 KiB key-page change showed little installed-R benefit at 10M and mixed paired timings. Root selected the R-only lifetime cleanup with the existing 16 KiB×512 cache. This review supports that bounded decision; it does not establish that 8 KiB pages are generally slower. The 8 KiB candidate passed correctness qualification and remains an archived, unselected experiment.

Timed role A loads qualified R-only B (16 KiB pages); timed role B loads combined C (8 KiB pages). Both have 8 MiB key-page payload. Values below use **100 × (combined_C − R_only_B) / R_only_B**: negative time/CPU is faster; positive RSS is higher. Percent summaries are medians of paired changes, not ratios of separate medians.

| Rows | Pair/order | R-only B op+closingGC s | Combined C s | Time change | Whole CPU change | Peak RSS change |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1M | 1/AB | 14.055 | 13.945 | -0.783% | -0.508% | +0.0327% |
| 1M | 2/BA | 14.075 | 14.014 | -0.433% | -2.250% | +0.1706% |
| 1M | 3/AB | 13.886 | 13.775 | -0.799% | -1.050% | +0.0716% |
| 10M | 1/AB | 131.622 | 131.417 | -0.156% | -0.972% | +0.0095% |
| 10M | 2/BA | 133.710 | 130.223 | -2.608% | -3.554% | -0.1279% |
| 10M | 3/AB | 130.112 | 132.709 | +1.996% | +3.088% | +0.0619% |

| Median paired change | Op+closingGC | Whole worker wall | Whole worker CPU | Peak RSS |
| --- | ---: | ---: | ---: | ---: |
| 1M | -0.7826% | -0.6813% | -1.0505% | +0.0716% |
| 10M | -0.1557% | -0.0150% | -0.9716% | +0.0095% |

All three 1M time pairs improved modestly. At 10M the changes were −0.156%, −2.608% and +1.996%; median operation-plus-closing improved only 0.156% and whole-worker wall only 0.015%. Median paired RSS increased at both scales: +753,664 B at 1M and +102,400 B at 10M. This demonstrates no memory improvement. The same nominal 8 MiB payload does not bound allocator/map/list overhead.

The prior, separate original 8b versus R-only trial at 10M showed 3.63% lower operation-plus-closing time and 4.68% lower whole-worker CPU, with all three pairs faster, while RSS increased 0.40%. That independently qualified R cleanup is the selected change. Its effects must not be added to or multiplied by this incremental experiment or earlier native/Darwin gains. The standalone 8 KiB native result (3.47%/3.78% faster matched/shuffled) remains valid for its own harness; it did not predict a comparable installed-R gain here.

## Qualification and measurement scope

One fresh candidate installation reproduced all 269 Linux cases / 3,077 assertions, zero skips/failures/errors and both existing warnings. Qualified R-only B was reauthenticated and reused. Twelve fresh timed workers ran three serial AB/BA/AB pairs per scale, with identical 100,000-row / 8 MiB batches, one-thread settings, references, options and unchanged worker code. A caller full closing collection is included in the reported operation total. Linux wait4 covers complete worker wall/CPU/RSS; output, RDS and oracle checks occurred outside timed workers.

Cloud verification covered 112 full stable-file comparisons (8 files × 7 editions × 2 scales), stable RDS fields, both raw-source copies, complete independent 19-field oracles on all 1M and 10M C pair-01 rows, transitive full-byte equality for the remaining outputs, namespace/dependency paths, physical installed trees, all 17 source files and 53 protected retained identities/hashes. The 10M gold is the retained qualified R-only pair-01-B edition; the 18 retired historical large files were not required.

The guard exited 0 at 2026-11-22T05:53:59.909259Z after 1,826.160 s with clean unsignaled closure. All 88 metadata files (905,074 B) returned, no absences. All six C workers were below 8 GiB (maximum 1,074,593,792 B). The separate full-suite/group peak was 5,998,723,072 B under its 24 GB cap. Final cumulative growth was 64,063,808,512 B under 80 GB, retaining the original pre-pack baseline and 186 GiB workspace floor.

The panel repeats 437,458 chr22 SNP loci with synthetic effects/frequencies, format 2, matched order. Three pairs provide no formal confidence interval or corpus-wide estimate. OS page cache was uncontrolled and scientific files were read for verification. Real-study/whole-genome, format 3 throughput, indel/multiallelic throughput, physical I/O, Windows, legacy format_sumstats and the historical cold/warm/40M gates remain unmeasured by this trial. Correctness qualification does not establish performance generality.

## Source, independent review and retention

- Experimental candidate: `8bc67611e96c588db79e853271951e2611056756`, tree `ab4ac6901713aceeabfc9d9e7fbd36db6e6a3bca`, direct parent `8b85aa5ffc45ad7632cfc5319e85dd5eb3089be6`.
- Selected R file SHA-256: `32f23a6779483741dc92642ccd0c15b859d5d485f25316934ae8a3159885e94a`.
- Retained 16 KiB header SHA-256: `5881233fbad917d4475f4801497c90a720428c41ba5afb8a7d38c234da60f155`; unselected 8 KiB header: `aff90b274a89f0063a6376011ddf585208c88ba954cb39c569cd0c105989a8f4`.
- Candidate source archive: 12,892,160 B, SHA-256 `da9ca1e47dbae6bdd4984d1f597a3e8a8b09a4f04887f28613f917d4e4a81363`; payload manifest `77c918632a63d1b19fff708ca8dcf025d1becb0c0fe914430b81bd98158177de`.
- R-only installed B tree: `d732c62becffbe538cf7250842bd73ee87c038225c4ee20cbf5cd79a381798e5`.

`review-analysis.py` independently reconciled all 47 metrics at each scale and 282 paired field comparisons using Decimal inputs, all medians, 88 returned pins, raw suite outcomes, gold/oracle/RDS/namespace receipts and all 24 operator / 17 payload postimages. `REVIEW.json` records the pass with no Critical/Important findings. `analysis.json` retains full precision and every GC component, including undefined zero-denominator percentages. `SOURCE-PINS.json` binds all local inputs/source/evidence and the 53 remotely protected hashes. No R/compiler/native/scientific or cloud call was made for this review; no repository integration or publication was performed.

The experimental source/archive and reviews were already independently retained (26 files / 13,094,464 B, zero differences) at `r2:smb-data-prod-scratch/gwas-harvest/mss-linux/2026-11-22/combined-installed-R/`. That source-retention checkpoint predates execution; runtime/result retention is a separate Root step. The unselected experiment and its measurements are preserved.
