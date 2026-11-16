# Scope correction for the preserved preliminary FAI evidence

The original `fai-slice-evidence.json` is preserved byte-for-byte. Its
`index_builds: 0` field means zero count-process index outputs and zero full
production index builds. It does not describe every test process.

The four preliminary FAI logs invoked the actual writer on76 tiny test fixtures:
12 initial exact-field missing-feature RED cases,40 writer calls in the41-case
GREEN suite,12 revised baseline RED cases, and12 final revised GREEN cases.
These tiny writer-oracle outputs are expected verification artifacts. They do not
admit a full source or production index. The failed private Fasta.contigs access
attempt made no production header/API change.

The new coordinator's final ordinary suite additionally creates42 tiny writer
oracles (40 original/FAI comparisons and2 comparisons against actual tiny producer
publications). Its strict sanitizer replay creates40 more. Earlier exploratory
new-coordinator reruns are logged separately and are not included in these final
per-run totals. Every count process remains read-only and creates no index output;
full production index builds and full-reference counts remain zero.
