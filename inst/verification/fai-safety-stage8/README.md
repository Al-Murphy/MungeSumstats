# Native FAI numeric safety and current-source qualification

This closes the inherited numeric-index follow-up recorded by the relocation
review. FAI decimal fields are parsed completely with checked uint64 arithmetic.
Malformed names/rows, duplicates, zero/impossible widths and last-base addresses
outside the opened FASTA are refused before native line arithmetic. Query ranges
use subtraction to avoid unsigned addition wrap. Valid leading-zero decimals,
UTF8 names, LF/CRLF/other consistent terminators, final lines without a newline
and the native writer's empty-contig representation remain supported.

This validates numeric/layout/address safety. It does not authenticate a newly
operator-rewritten manifest/index against every FASTA byte or claim to reconstruct
the whole sequence layout. Original trusted bundle checksums still decide identity.
The successful normalization code following reference handling is byte-identical
to the reviewed baseline; no R/API/generated wrapper/scientific version changes
were made. Reference formats2/3 and optional relocation retain their interfaces.

Production candidate **03d11e92b417add0183dea4fb77716f5502d5735**, core SHA256
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`.
Independent review is **GREEN**. Its11,271 constructor cases matched a separate
big-integer oracle (741 accepted/10,530 refused);35,535 exact sequence queries,
including160 valid physical layout controls, passed ASan/UBSan without findings.
Public R probes, saved installed evidence and full source/DLL integrity were checked.
Earlier baseline mismatches, test-only assertion corrections and unsupported host
diagnostic-wrapper failures remain retained rather than counted as successful tests.

[Native CI43858052096](https://github.com/Al-Murphy/MungeSumstats/actions/runs/43858052096)
passed Windows, macOS and Linux at5724713 (CI branch trigger only beyond03d11e9).
Each OS ran82blocks/1,727SourceCpp assertions and42blocks/814assertions through
the actual generated registered bridge. Unix-only standalone native fixtures are
additional; that shell step is not claimed on Windows. Full installed package:
**257blocks/2,570assertions**, zero failures/errors, two known warnings and two
missing-reference skips. Exact installed DLL:
`6f23605e06ac0b35fd40268b08190fd55e9a09e2fccded2f8cd7b629a3cd49eb`.
All868 archived source files and the DLL were rechecked afterward. The installed
suite passed its8GiB/160GiB-free/4GiB-output/900second guards in179.523seconds;
peak owned-tree RSS6,055,772,160B and physical5,831,824,760B. This suite peak is
not a whole-file processing memory measurement.

An append-only profile ran the unchanged independently reviewed whole-file harness
against this source/DLL. One R/data.table thread, fresh application handle/process,
shared Apple M5 Pro host, uncontrolled OS page cache, reference format2, no projection.
One million rows repeat437,458 unique real nonpalindromic biallelic chr22 SNPs with
synthetic exact effects/frequencies; this is not one million unique variants or
real-study coverage. Indel/multiallelic correctness is covered by the separately
reviewed scientific matrices, not timed by this workload.

| Matched1M case | Processing | External wall | External peak RSS |
| --- | ---: | ---: | ---: |
| Earlier preservation implementation4 |20.211s|24.900s|1,617,166,336B|
| Scientific7 plus this FAI safety change |24.444s|29.010s|1,629,683,712B|

Every row/all19 normalized fields and both original/source copies passed the
independent oracle; all1M normalized rows are byte-identical across cases, with
zero unexplained deltas. Match/swap counts are562,542/437,458. Both process guards
and the separate verifier passed; minimum observed free disk466,825,166,848B.
The current stricter scientific backend is about20.9% slower than implementation4
on this workload, with similar memory. This comparison is not against legacy
`format_sumstats`. Earlier legacy-memory/reference-preparation gains remain in the
prior archives with their original source pins. Earlier10M/40M timings are not
relabeled for this source; no larger-stage benchmark was executed here.

`archive-map.json` binds raw and archived byte counts/SHA256; logs and decoded CI
JSON are losslessly gzip compressed. The unchanged harness/seal, exact new profile,
root1M-only approval, source-bound installed receipt and per-case raw timing/guard
logs are included. `omissions.json` pins generated binaries/editions/reproducible
case trees; generating programs and original complete evidence indexes remain.
Large inputs/references are external immutable pins, not vendored. No actual cloud
credentials, account inventories or active full-reference files were touched.
