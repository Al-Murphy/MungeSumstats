# Scientific MSS7: complete whole-file streaming benchmarks

Frozen scientific source: `5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e`; DLL `c1c0a8b29479c33d8b4cfa4646e46e39c4da436686f4ca2da12a402251003d50`. Reference format2; one R/data.table thread; projection not requested. The optional typed-format3 prototype and its measurements have separate source/DLL pins.

**All stages completed and passed independent, exhaustive output verification.** Each run retains two byte-identical source copies and checks all19 normalized fields of every row. Across the matched1M pair and scientific7 10M/40M runs,52million normalized rows were examined, with zero unexplained deltas.

|Rows|Implementation|Processing (s)|Whole process (s)|Peak RSS (GiB)|Peak physical footprint (GiB)|
|---:|---|---:|---:|---:|---:|
|1,000,000|mss4|20.341|24.96|1.520|1.308|
|1,000,000|mss7|24.192|28.89|1.499|1.251|
|10,000,000|mss7|256.408|261.07|1.808|1.350|
|40,000,000|mss7|974.111|979.09|2.131|1.390|

The current matched million-row pair measures **18.93% processing overhead** for scientific7 versus frozen4. Peak RSS is similar. The added ambiguity/context checks buy scientific correctness; this workload shows their CPU cost rather than an ordinary SNP throughput gain. The40M run processes in974.111seconds with2.131GiB external peak RSS and completes within the8GiB process-tree target. Its sampled output/temporary high-water is9,413,637,648bytes; the external shared reference is separately pinned.

## Workload and correctness boundaries

The inputs repeat437,458 unique real chr22 SNP records with synthetic beta/EAF, alternating source effect orientation in full-panel cycles. The40M input is **40million rows, not40million unique variants**. The reference was confirmed to contain437,458 unique biallelic non-palindromic loci/primary IDs. The independent oracle refuses duplicate, multiallelic, complex or palindromic references for this deliberately narrow performance workload. Indel/multiallelic identity, direction, frequency and uncertainty behavior is covered by the separate scientific7 mathematical/public-R reviews.

The matched1M outputs are byte-identical across all19 normalized fields. The40M outcome counts are20,123,068 match and19,876,932 swap; effects and EAF are independently exact for all40M rows. This does not estimate complete real-study normalization rates or production-reference coverage.

## Measurement and reproducibility

Each case uses a fresh R process and verified reference handle. Package setup, reference verification/open, processing, external process time and correctness-verifier time are distinct receipt fields. OS cache is uncontrolled and never flushed; the shared host also runs independent reference preparation and review workloads. These are single measured trials, not confidence intervals or production service latency promises.

The frozen v2 harness has independent GREEN review,14 adversarial verifier selfchecks, exact source/library/reference/input pins and an explicit coordinator science approval. Version1 was not executed: its approval-text/row-width/reference-uniqueness review findings and RED receipts remain preserved. A separate test-only path-alias assertion was corrected before sealing v2.

All owned commands are guarded at8GiB tree RSS/physical accounting,160GiB free-space floor, finite deadlines and explicit owned-output growth caps. All benchmark and verifier guards succeeded; no censor/OOM occurred. Sampling can miss sub-poll transients; external time separately reports process peaks. The full40M verifier itself took294.83seconds and used about1.47GiB sampled RSS; that cost is not hidden inside harmonization processing.

Reproduce from the frozen harness by regenerating a fresh pinned profile from the same source/reference/input bytes, obtaining an exact source/DLL scientific approval, then running matched1M followed by authorized larger stages. Original absolute paths/timestamps are forensic receipts; do not edit a sealed profile to relocate it. The fixture generation/acquisition recipe and historical `format_sumstats` trials retain their original separate pins.

Full GRCh37/GRCh38 references, complete real-source MSS harmonization and actual full1000 execution remain independent pending gates. No full-corpus capacity, legacy40M success ratio, universal harmonization or typed3 whole-file performance claim is made by this format2 report.
