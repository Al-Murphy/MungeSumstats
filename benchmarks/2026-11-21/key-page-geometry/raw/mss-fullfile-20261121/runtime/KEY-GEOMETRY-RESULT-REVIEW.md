# Key page geometry independent result review

**The native mechanism gate passes. Critical 0, Important 0, Minor 0.** The 16 KiB trial merits a guarded unprofiled public-API comparison. It preserved the native results and reduced nominal requested key bytes in both orders at the same page payload budget. More key-cache misses increased application stream seek/read calls, so the elapsed-time effect remains open. The earlier minor finding on numeric diagnostic context stays deferred to promotion work. [Source review](KEY-GEOMETRY-SOURCE-REVIEW.md)

## Run and artifact verification

- Root observed outer session `31893` exit `0` in tool chunk `2dbd6f`. I verified the observation SHA-256 `64ecf5be7a7827859ec5e702f4e38e17cedcb5032de447fa4c408e96853022cd` and each of its eleven guard receipt pins. Every receipt has a planned, running and completed record with raw exit `0`, confirmed closure, a reaped leader, empty stderr, and no trigger, signal or cleanup error. The observation records all eleven groups and sessions absent. [Root observation](/private/tmp/mss-key-geometry-experiment-20261121/root-observation-01.json)
- Both guarded compile plans used Apple clang 21.0.0, the macOS 27.2 SDK, `-std=c++17 -O2 -arch arm64 -DMSS_LOOKUP_COUNTERS`, the same driver and assertions enabled. The plans differ only in role-specific include and executable paths. I verified the staged header pins and saw only the reviewed change to the keys initializer. The executable SHA-256 values are `1cab768d53d4b44f1fd4d6c11bbb1fad4093e155772a8f93002102fd9700ab46` and `aa6127ba1cb59b61f3827d1a9cd75c0b425d338428beadd7ed7a7a1898fbd37d` for baseline and trial. The product checkout reports `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b` with no tracked changes. [Admission](/private/tmp/mss-key-geometry-experiment-20261121/run-01/admission.json) · [Trial header](/private/tmp/mss-key-geometry-experiment-20261121/run-01/trial/src/preservation_core.h:502)
- All four matched plans point to the same pinned format-2 index, FASTA and input, with distinct role and order outputs. The seven recorded file identities agree before and after and match the admission record. The receipt pin is `873e8d12793460bd351cf2521d1cfc17db7cf2b7d3738ea2f6646d0eb16ffae0`. I checked the receipt and did not rehash fixture payloads. [Input window](/private/tmp/mss-key-geometry-experiment-20261121/run-01/input-window.json)

## Native parity and selection

The legacy parity stdout bytes match across roles for formats 2 and 3, with SHA-256 `57da31b533f1b9996ab5a76edc0f4930a4c0cc1abe7b622c4a982cf624f6ba4e`. Both pressure invocations passed with assertions enabled. Both formats and roles report zero allocation calls for the 1,024 warm-negative queries and complete the 3,000-distinct-query shuffled workload. [Baseline parity](/private/tmp/mss-key-geometry-experiment-20261121/run-01/parity-baseline.stdout) · [Baseline pressure](/private/tmp/mss-key-geometry-experiment-20261121/run-01/pressure-baseline.stdout) · [Trial pressure](/private/tmp/mss-key-geometry-experiment-20261121/run-01/pressure-trial.stdout)

I streamed 18,324,100 bytes across the four runs' outcomes and order files. The output and order bytes match line for line between roles within each order. The 21-column schema contains the source ordinal, nine raw tokens and all eleven native fields, including both serialized 64-bit float patterns. The input ordinals equal the prescribed floor formula with endpoints 1 and 437,458. The shuffle is a different permutation of exactly those ordinals. Joining by ordinal gives identical records across orders. Each run contains 32,768 distinct raw and 32,768 distinct canonical tuples, and every selected event resolves as `match`. The real panel therefore exercises successful ordinary matching. Other orientations and malformed inputs are represented by the legacy characterization fixtures. [Baseline matched input](/private/tmp/mss-key-geometry-experiment-20261121/run-01/baseline/matched-input/outcomes.tsv) · [Baseline matched shuffle](/private/tmp/mss-key-geometry-experiment-20261121/run-01/baseline/matched-shuffled/outcomes.tsv)

| Order | Outcomes SHA-256, both roles | Ordinal order SHA-256, both roles |
|---|---|---|
| Input | `18084704d5328e1e68b6ccf41b7c284a206110386866f6dab989f2c832a543e4` | `d79a024309196f45682001d844a8eef2cbbbeb61fc06f43dc350540df8fd7d1a` |
| Shuffled | `14cd319058bdc304eff43321fe9343e07fbb376575be4b9e6b32517c22d85bb8` | `7d117f0390a21530979d0482531bd4e817e205fc9a2546199820d5da61f1082a` |

## Cache mechanism

| Order | Key page reads, baseline → trial | Read change | Nominal requested key bytes, baseline → trial | Byte change |
|---|---:|---:|---:|---:|
| Input | 4,746 → 7,032 | +48.167% | 311,033,856 → 115,212,288 | −62.958% |
| Shuffled | 73,336 → 78,265 | +6.721% | 4,806,148,096 → 1,282,293,760 | −73.320% |

The `counters.tsv` values agree with `result.json` (SHA-256 `26150d11b7078b18f412ff8ff35ac05fbc9d3a64d10ad379d8e881ea6593e0bc`). Key page requests increased from 932,966 to 934,690 in each order (+0.185%). All other native counters matched by role within each order. Both orders have 735,785 probe entries, 1,269,761 pivot hits, zero rank hits, 932,392 offset requests, 8,191 pivot slots, 4,096 rank slots and 983,056 rank bytes. Offset reads were 164 in input order and 20,712 in shuffled order. [Run result](/private/tmp/mss-key-geometry-experiment-20261121/run-01/result.json)

`BlockFile` counts a page request at cache entry and a page read on a miss after the stream seek and read. Multiplying miss counts by the configured page widths of 65,536 and 16,384 bytes gives the table's requested byte totals. The stream can return fewer bytes on its final short page. Filesystem caching and buffering can independently change physical I/O. The trial trades smaller reads for more application stream seek/read pairs, especially in input order. The same probe, pivot, rank and offset counts isolate the observed change to key page handling for this panel. [Page cache implementation](/private/tmp/mss-key-geometry-experiment-20261121/run-01/baseline/src/preservation_core.h:36)

The staged headers use 64 KiB × 128 and 16 KiB × 512, both 8 MiB of key page payload per handle. Offsets and records budget 8 MiB each and FASTA 4 MiB, for a total of 28 MiB in both roles. Page container overhead may vary. Guard RSS had only two to four samples per matched phase and describes the standalone processes. Public-API RSS and paired unprofiled elapsed time require separate measurement. [Default cache geometry](/private/tmp/mss-key-geometry-experiment-20261121/run-01/baseline/src/preservation_core.h:69) · [Trial initializer](/private/tmp/mss-key-geometry-experiment-20261121/run-01/trial/src/preservation_core.h:502)

## Next gate

Run a guarded, unprofiled baseline/trial comparison at the normal matched public API. Compare the complete 19-field oracle, measure elapsed time and peak RSS across paired runs with controlled ordering and cache conditions, and require a repeatable benefit before promotion. Keep format-2/3 and key-page boundary characterization and native Windows qualification in the promotion scope. This review ran only bounded artifact checks. [Experiment plan](KEY-GEOMETRY-PLAN.md)
