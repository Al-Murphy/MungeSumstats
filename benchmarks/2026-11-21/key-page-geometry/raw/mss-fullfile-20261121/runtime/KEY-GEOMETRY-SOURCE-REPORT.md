# Key page geometry source handoff

The scratch driver adds `lookup_parity matched INDEX FASTA INPUT OUTPUT ORDER`. The trial patch changes one initializer in `Reference` from `keys(dir+"/keys")` to `keys(dir+"/keys",16384,512)`. The product checkout at `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b` remains at the same revision, with no tracked changes. Compiler and guarded execution results are pending.

The candidate compares 64 KiB × 128 with 16 KiB × 512, both 8 MiB of key payload per handle. The patch leaves the default `BlockFile` configuration, offsets and records at 8 MiB each, and FASTA at 4 MiB, for 28 MiB total page payload. The trial can hold four times as many key-cache entries, so cache container overhead can differ. The experiment targets the lower-bound key-line stream-read path, which accounted for 514 of 3,446 main-thread observations in the current sample. Page counts will measure application cache reads. Multiplying key reads by the respective page width derives nominal requested bytes, which can exceed returned bytes for a short final page. Physical I/O, elapsed time, and RSS require separate observations. [Profile diagnosis](CURRENT-MSS-PROFILE-DIAGNOSIS.md)

## Artifacts and pins

| Artifact | SHA-256 |
|---|---|
| [Committed original driver](/private/tmp/mss-lazy-variantannotation-20261121/tests/standalone/lookup_parity.cpp) | `85dd1ac49ae84089504708c212ef57463cafecab52943f9cd5717a0c595a38c1` |
| [Committed baseline header](/private/tmp/mss-lazy-variantannotation-20261121/src/preservation_core.h) | `42cbd86be7e8e85dd67b11713db3f3483997be96f6ace69b54d655d0ca7decf5` |
| [Scratch matched driver](/private/tmp/mss-key-geometry-experiment-20261121/source/lookup_parity.cpp) | `505b0e8c72a517e7db4b48fb49126c0bfcbeeb6254dcd19c8090db57dd20fec9` |
| [One-line trial patch](/private/tmp/mss-key-geometry-experiment-20261121/source/trial-header.patch) | `cdfc635afd0c1b7f654301801c12d0f5bbab4393da830b9b0e259ebef81c96cf` |
| Expected trial header from the in-memory initializer substitution | `5881233fbad917d4475f4801497c90a720428c41ba5afb8a7d38c234da60f155` |
| [Full driver and header review diff](KEY-GEOMETRY-SOURCE.diff) | `27853e5bee84e0bbee8a1464ac801ceb5344a8edcb76e195f767ffc03c9ff097` |

The invocation uses the real [format-2 matched index](/private/tmp/mss-fullfile-20261121/fullcall-1m-01/matched-index), [chr22 FASTA](/private/tmp/mss-fullfile-20261121/asset-restore-01/GRCh37.ensembl104.chr22.fa), and [matched input](/private/tmp/gwas-harvest-stage2/assets/mss/matched-1M.tsv). The input SHA-256 supplied by the [experiment brief](KEY-GEOMETRY-BRIEF.md) is `d110828a8b17bf88309db042e7073faa267a0b63a3550d14a87f7761074c91f2`. Input authentication for the execution gate remains with the guarded run owner.

## Matched mode contract

The driver reads the first 437,458 data rows with the exact nine-column header and reuses `mss::split` for this fixture's literal-tab fields. It checks each row's field count and framing, a fully parsed positive `BP`, and fully parsed finite `BETA`, `FRQ`, `SE` and `P`. It retains the original line for each selection, preserving the lexical tokens including decimal spellings. The selection uses 1-based source ordinals `1 + floor(i * (437458 - 1) / (32768 - 1))` for `i = 0..32767`.

Selected raw event identity is `(CHR,BP,A1,A2)` using source tokens. The driver maps `CHR,BP,SNP,A1,A2,BETA,FRQ` into the ordinary `Input` with `A1` as other allele, `A2` as effect allele, and frequency type `EAF`. It opens one fresh `Reference` per invocation and calls `harmonise` once per selected row with candidate cap `10000` and the default disabled source contract. Each returned row must have status `match`, `swap`, `complement`, or `swapped_complement` and a distinct canonical `(chr,pos,ref,alt)` tuple.

`input` visits the selected rows in increasing source ordinal order. `shuffled` uses `std::shuffle` on the selected indices with `std::mt19937(257)`. Both builds use the same driver and toolchain for each order. The permutation applies to that build environment. Stdout records the order, seed, selection formula and both event identities.

The output directory must be fresh. Files use literal tabs and LF in binary mode. The start of each outcome record is `source_ordinal` followed by `SNP`, `CHR`, `BP`, `A1`, `A2`, `BETA`, `FRQ`, `SE`, `P`. The native columns follow in order as `status`, `chr`, `pos`, `ref`, `alt`, `id`, `effect_bits`, `eaf_bits`, `candidate_count`, `swapped`, `complemented`. The existing `bits(double)` emits both floating-point patterns as unsigned decimal integers and the two flags emit as `0` or `1`. `order.tsv` records one source ordinal per query in execution order. `counters.tsv` uses `metric` and `value` columns and records numeric values for `selected_count`, `distinct_raw_count`, `distinct_canonical_count`, then `rank_hits`, `pivot_hits`, `probe_entries`, `key_requests`, `offset_requests`, `key_reads`, `offset_reads`, `pivot_slots`, `rank_slots`, `rank_bytes` from the actual pass. It is written after the outcomes and order streams close successfully. A failed invocation can leave partial output in its scratch directory. The successful exit and the complete outputs are the comparison inputs.

## Source checks

- Read the driver and header bytes from the named Git commit and confirmed equality with the checkout files. Removing the matched helpers, matched `main` dispatch and two added includes from the scratch driver reproduced the original bytes, covering the legacy `parity` and `pressure` bodies and dispatch. The original `NDEBUG` guard remains, and the matched dispatch requires `MSS_LOOKUP_COUNTERS`.
- Both the trial patch and the full review diff passed `git apply --check --whitespace=error-all` against the checkout. The full diff reports 149 added driver lines and one header line replaced. Independently checked that the header hunk adds only the two key-cache arguments and computed the expected trial header pin from an in-memory substitution.
- Inspected the matched source against the native `Input`, `Result`, `LookupCounters`, and ordinary resolution statuses. A source check verified the declared and emitted counter names, input mapping, event tuple expressions, native outcome column order, and the two `bits` calls. The first checker pattern also picked up `source_ordinal` from the outcomes header. Scoping that pattern to the counter stream gave a clean rerun with the source pins unchanged.
- Checked the actual fixture header bytes and evaluated the ordinal formula independently. It yielded 32,768 strictly increasing positions from 1 through 437,458, with gaps of 13 or 14. Resolution, raw and canonical distinctness, native parity, and counter values will be evaluated against the pinned panel during guarded execution.

The source checks did not invoke a compiler, R, or the native driver. Whole public-API and 19-field oracle parity, unprofiled timing and RSS, new page-boundary coverage, and native Windows validation remain conditional promotion gates.
