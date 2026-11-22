# R batch-lifetime cleanup

This revision releases decoded file-batch intermediates after normalization and releases each normalized table after writing and counting it. It retains the existing 16 KiB × 512 key cache. The selected R file has SHA-256 `32f23a6779483741dc92642ccd0c15b859d5d485f25316934ae8a3159885e94a`; the native header remains `5881233fbad917d4475f4801497c90a720428c41ba5afb8a7d38c234da60f155`.

The exact R change passed the complete installed Linux suite (269 cases, 3,077 assertions, no skips/failures/errors, two existing warnings), complete output/RDS parity, and independent 19-field row oracles at both scales. Three fresh serial AB/BA/AB pairs compared against installed 8b with identical one-thread, 100,000-row / 8 MiB settings.

| Rows | Median paired processing + closing-GC reduction | Whole-worker CPU reduction | Peak RSS change |
| --- | ---: | ---: | ---: |
| 1M | 2.02% | 1.57% | 0.45% lower |
| 10M | 3.63% | 4.68% | 0.40% higher |

All three processing-plus-closing pairs improved at both scales. This is a modest CPU cleanup, with no sustained RSS reduction demonstrated. The closing collection is included. GC timings are durations, not collection counts or heap occupancy.

The 1M producer originally exited 1 after its six timing workers and full suite passed: its verifier omitted a defaults merge. A separate verification-only completion authenticated the untouched outputs, source and complete row oracle without rerunning timing. Both historical outcomes remain preserved. The 10M trial completed normally. The original reports below retain their historical status wording.

A later independently reviewed installed-R experiment tested 8 KiB pages on top of this R cleanup. Its incremental processing-plus-closing reductions were only 0.78% at 1M and 0.16% at 10M, with mixed 10M pairs and flat RSS. It was not selected; see [the complete experiment](../key8k-installed/RESULTS.md). These experiments and older native/Darwin results must not be added or multiplied.

These panels repeat 437,458 real chr22 biallelic, non-palindromic SNP loci with synthetic effects/frequencies. They do not establish whole-genome, indel/multiallelic or real-study throughput, Windows or legacy `format_sumstats()` speedups. Full historical cold/warm/40M and broader source calibration remain open. No MSS scientific decision, public interface, dependency, reference format or compiler flag changes in this cleanup.

Original source/evidence and failed/completed verification receipts were independently retained under the operator-configured R2 prefix `gwas-harvest/mss-linux/2026-11-22/`. The later experiment's complete 123-file source/runtime/report bundle is under `combined-installed-R/evidence/`, with zero differences on independent download check. Its source archive was retained separately. Content hashes bind the evidence; storage account names are not part of package behavior.
