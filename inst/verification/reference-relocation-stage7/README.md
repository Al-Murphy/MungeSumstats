# Verified FASTA relocation and native portability

`open_reference(reference, fasta = runtime_path)` can explicitly bind a transferred
reference to a byte-compatible FASTA at its new location. It rereads the original
manifest, verifies its original FASTA size/MD5 and every immutable reference artifact
including source and FAI, and opens a fresh native handle. It never guesses a path,
rewrites the manifest, skips verification or reuses a serialized native pointer.
Physical binding receipts are separate from scientific checksums/versions. Default
calls and default reference format2 retain their prior behavior; format3 is opt-in.
Harvest separately pins the transferred bundle with SHA256.

Production source is the independently reviewed `8334f1c` R-only change. The final
native CI head is `50912b79d93ac7ee14127cb3a83fb4c16a5060d5`; its R/src bytes are
identical to8334. Later commits adjust CI coverage and portable test fixtures only.
The installed DLL is
`b63e4e43a167ad5b4807bc6a36758e820eb7b0c70ae3e77124e6ab909d9e4841`.

[Native CI43848950725](https://github.com/Al-Murphy/MungeSumstats/actions/runs/43848950725)
passed on **Windows, macOS and Linux**. Each OS ran80blocks/1,427assertions through
SourceCpp and40blocks/514assertions through the actual registered generated bridge.
The registered path compiles `preservation.cpp` plus `RcppExports.cpp`, resolves
all seven registered entry points and sources the real generated R wrappers.
These ordinary CPU/R tests require neither WSL/containers nor external bcftools.
The additional standalone Unix C++ gold-fixture step is still skipped on Windows;
the two R/native test steps run there.

The final full installed suite passed **255blocks/2,270assertions**, zero failures
or errors, two known upstream warnings and two missing-reference skips. Guarded
wall time153.553seconds; peak tree RSS6,099,501,056B, physical5,869,262,440B;
minimum observed free disk474,044,379,136B. The8GiB/160GiB/4GiB-output/900second
guards passed, and234 tracked source hashes plus the DLL were rechecked afterward.
The earlier255block/2,262assertion run is preserved separately.

Earlier CI attempts are retained with their exact heads and decoded logs. Windows
refused a timestamp change while a native FASTA stream was open: the old fixture
ignored `Sys.setFileTime(FALSE)` and incorrectly expected mutation detection.
The corrected fixture asserts the actual write outcome, then changes exactly one
binary byte of an unheld source artifact and checks a fresh binding's refusal.
An intermediate text-mode LF append became CRLF on Windows; it was replaced by
binary append. Offline directory moves now explicitly collect unreachable native
streams left by completed streaming calls before moving files. All source/normalized
parity, immutable-manifest and real-mutation controls remain exercised; no platform
skip was added to hide these failures.

Original independent review and the subsequent coordinator fixture/worker review
are **GREEN**. The latter independently reran13actual worker cases with each library,
both local native paths and an8block/192assertion forced timestamp-denial control;
it independently verified the saved native CI and installed-suite/source pins.
Harvest worker adoption and its combined integration are separate from the R change.
One inherited native-parser follow-up
remains: an operator-rewritten, newly trusted manifest can bind a syntactically
five-column FAI with invalid zero line widths. Ordinary corruption against the
unchanged manifest is refused. This is outside this R-only relocation change and
is not a bypass of the unchanged trusted bundle hashes.

The scientific7/format2 whole-file benchmarks and typed3 reference stress results
remain in `../reference-format3-stage7/`. They retain their original source/DLL and
workload pins; they are not relabeled as relocation-source throughput measurements.
This relocation suite does not establish full production-reference completion,
full-study normalization coverage, or public cloud readiness.

`archive-map.json` maps original paths to archived files, with original/archived
SHA256 and byte counts. Logs and decoded Actions JSON bundles are losslessly gzip
compressed. `omissions.json` pins generated RDS results omitted in favor of full
CSV/summary/log receipts; installed libraries and build binaries are not vendored.
Native receipts contain decoded API logs, not the original GitHub ZIP byte stream.
Original evidence indexes keep their original absolute paths; resolve them through
the archive map. No actual cloud credentials or private AWS inventories are included.
Unrelated hook-generated memory entries present in a reviewer index are explicitly
omitted without changing that original index; they are not verification artifacts.
