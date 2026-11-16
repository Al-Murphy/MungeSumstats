# Sequential FASTA offsets: measured MSS reference preparation

The portable C++ FASTA scanner counts physical bytes as it reads each line,
replacing two `tellg()` calls per line. It counts CR before stripping it and
counts LF only when `getline()` consumed a delimiter. Normalization, identifier
resolution, reference format and the legacy `format_sumstats()` API are unchanged.

## Measurements

All runs used the same complete derived GRCh38 FASTA: 3,339,671,278 bytes,
705 contigs, SHA256
`d88cc10556244b91b044294067ad166d9512bec6bbc104fb5c9f2d480af302fe`.
The index input was explicitly incomplete: the first 1,000,000 sites, without
complete merge-history expansion. Its exact 31,408,384 bytes have SHA256
`3136a1b6bb0dd7e6bba90410bacf47ac2e10eaab0d07e190a38a1be8ed066372`.

| Measurement | Frozen baseline | Candidate | Observed ratio |
| --- | ---: | ---: | ---: |
| Standalone FASTA scan | 55.630 s | 2.525 s | 22.03x |
| Matched whole reference preparation | 140.709 s | 30.273 s | 4.65x |
| Whole preparation process-tree peak RSS | 1,050,902,528 B | 1,066,106,880 B | No memory improvement claimed |

The matched pair used the same `Rscript --vanilla` program, Rcpp preload, input,
batch size, reference version and flags, with separate frozen installed package
libraries. Baseline and candidate ran sequentially. The host was shared with the
ongoing bioDB reference build; OS cache was not controlled. These are observed
single-pair measurements, not isolated cold-disk estimates. They measure reference
setup, not full-study harmonization or a legacy-pipeline speedup.

All seven scientific payloads are byte-identical: `fasta.fai`, `keys`, `offsets`,
`records`, `rejected.tsv`, `source.tsv`, and `core-version`. They total
232,818,142 bytes. `manifest.rds` has operational path/time differences and is not
claimed byte-identical. The external FASTA is additional. Six unsupported source
alleles remain audited in `rejected.tsv`; their site context remains indexed.
The index-size sample does not establish full-genome or complete-history costs.

## Verification

* Installed complete suite: **204 test blocks, 1,268 passed assertions, zero
  failures/errors, two known warnings, two missing-reference-package skips**.
  Guarded wall time 141.390 s, peak process-tree RSS 6,389,121,024 bytes. The two
  skips are dbSNP144 GRCh38 and dbSNP155 GRCh37 packages. The warnings are the
  existing S4Vectors import replacement and rsID coercion cases, recorded in the
  full log. This is separate from the reference-preparation memory figure.
* The new persisted tests cover 24 binary LF/CRLF/terminated/unterminated/line-width
  layouts, public per-base resolution and source preservation, plus two single-line
  unterminated FASTAs. The first 144 assertions passed against the unchanged
  installed baseline before the scanner edit; all 146 pass in the candidate suite.
* Independent original/candidate C++ drivers checked 220 valid layouts and
  **173,355 sequence queries** against source-byte oracles. Twenty adversarial
  cases gave the same nine errors and eleven successes. The review is GREEN.
* The standalone R harness includes the new tests. Native Windows/macOS/Linux CI
  passed as recorded below.

Preserved setup diagnostics include an initial CRLF-header assumption failure
in the operator sample script, and a standalone invocation without the required
local Rcpp library path. Corrected configured runs passed. Namespace-only loading
of the old installed package emitted a stack-imbalance warning before native
calls; the matched vanilla/Rcpp-preloaded controls and installed suite did not.
This change does not claim to fix that separate loading behavior.

FASTA inputs still require fixed-width sequence lines. Mixed LF/CRLF inside one
contig can give inaccurate indexes in both versions; independent review confirmed
this pre-existing unsupported-layout limitation. Uniform LF/CRLF and different
endings between contigs were verified. The production derived FASTA has uniform
layout. No source/reference files were modified, and no full production MSS
reference or 1,000-input pilot was certified by this sample.

## Native CI receipt

[Run 43776925410](https://github.com/Al-Murphy/MungeSumstats/actions/runs/43776925410)
passed on native Windows, macOS and Ubuntu at source/docs commit
`a167b192137a31c684fbac1b4f82f6a4318c0670`. Each platform passed 29 R test blocks
and 425 assertions, including both new FASTA-layout blocks. The separate C++
gold fixtures passed on macOS and Ubuntu; that shell step is intentionally
not run on Windows, where the R interface itself compiled the portable C++ core.
See `native-ci.json`, `native-ci.log` and `native-ci-counts.json`.
