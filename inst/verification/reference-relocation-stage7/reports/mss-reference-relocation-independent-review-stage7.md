# Independent review: optional verified FASTA relocation

**Verdict: GREEN for the narrow R-only relocation candidate.** No new blocking
finding. Full installed-package validation and actual native Windows/macOS/Linux
CI remain coordinator gates; this review does not claim them completed.

Reviewed candidate `8334f1c9052ae266dc9df80a35bca07868af1423` against
`f1cd89cf5e0f24fa754bbe2f0e17c9d93c53d6a6`. The review ran from an independent
`git archive` export at
`results/mss-reference-relocation-independent-review-stage7/candidate`.
All **574** archived files matched the frozen commit after the runs; the author
checkout remained clean, `git diff --check` passed, and all native source and
generated native/R registration wrappers were unchanged from the base.

The supplied author report SHA-256 was independently confirmed as
`df288b3882d426bd33db1054dcdc4643726ba04e58d6f7a4c480662ad9431c99`;
its evidence index was
`9ac7a85ed0d27b9f46cbf308546b8783f19515a019d2e8a289a407ce1e4fd864`.
Every file named by that index and the indexed report was checked against its
size/hash. Those checks supplement the independently compiled and executed runs.

## Verification executed independently

| Run | Result |
|---|---|
| Fresh SourceCpp standalone suite | **80 blocks / 1,419 assertions**, exit 0 |
| Fresh generated registered `.Call` suite | **40 blocks / 506 assertions**, exit 0 |
| Additional independent probes | **8 blocks / 140 assertions**, no failures, errors, warnings or skips; exit 0 |
| Rd parse | Passed |
| Frozen source/artifact integrity | Passed |

Both native runners built the frozen source locally, rather than depending on an
author-built shared library. The registered runner compiled the actual
`preservation.cpp` and `RcppExports.cpp`, resolved the seven registered names and
loaded the real `R/RcppExports.R`. The review ran on macOS arm64 with R 4.5.3 and
local Rcpp/testthat; no Bioconductor package load or external bcftools was used.

The independent fixture includes SNPs, indels, a multiallelic source site, current
and historical IDs, a UTF-8 synonym, contig aliases, coordinate-only and ID-only
selectors, opposite measured-allele orientation, and exact source tokens. It uses
ordinary file copies under spaces and UTF-8 paths with a different FASTA basename.
No symlink privilege is required. Both format 2 and opt-in typed format 3 were
exercised; projection checks cover mixed 2→3 and 3→2.

The additional probes establish:

- `open_reference(reference, fasta=NULL)` preserves the original first argument
  and default. Existing handle reuse returns the same handle. Baseline R functions
  sourced independently from the base commit return identical complete default
  in-memory results. All nine ordinary edition files and eleven projection edition
  files are byte-exact; their complete manifests match after removing only output
  location and creation timestamp.
- An explicit replacement opens a different native pointer and uses the normalized
  actual physical path. It works after the original FASTA/index directory is absent
  when the caller explicitly supplies the new index path. An ordinary open without
  the required override refuses; no basename or environment inference is introduced.
- Rebinding a poisoned handle discards its altered in-memory checksums, sizes,
  files, FASTA path, build, contig aliases, signature and binding receipt. Rebinding
  a serialized handle creates a fresh native pointer; direct use of its dead pointer
  refuses. A real manifest checksum/version change after the first handle was opened
  is reread and refused on explicit rebinding.
- A same-size wrong FASTA refuses. Ordinary corruption of FAI, keys, offsets,
  records, source, contig-alias table and format marker refuses against the original
  manifest. Restoring the original bytes restores successful binding.
- Changes after verification and after native opening are both caught by the outer
  signature comparison. Runtime mtime changes invalidate reuse with `verify=FALSE`.
- Source-only and target-only projection relocation add only their applicable binding
  receipt. All pre-existing scientific receipt fields remain identical; complete
  normalized/source results remain identical. Mutation of either actual projection
  side refuses, including when harmonisation uses the original other side.
- Empty streamed input retains its explicit binding receipt and original reference
  checksums. Every unaffected original/transferred bundle hash, including RDS
  manifest bytes, remains unchanged.

Static review confirms physical binding fields are appended only to operational
table/edition/projection provenance. The existing scientific reference checksum
values, chain MD5, parser/normalization/projection versions and reference writer
logic remain unchanged. Ordinary edition `checksums.tsv` naturally changes for an
explicitly relocated edition because its operational metadata changes; that does
not alter the scientific reference checksums.

## Inherited FAI semantic validation follow-up

This is **not a relocation blocker or an immutable-bundle checksum bypass**.
The independent discovery probe deliberately replaced both a tiny fixture FAI and
the fixture's trusted manifest size/MD5. A four-column FAI correctly refuses in the
native parser. A five-column `1\t16\t3\t0\t0` FAI opens under both the candidate
and the base R API after the trusted manifest is rewritten to bless it. The native
constructor checks field count and numeric conversion but does not reject zero
line widths. No harmonisation was attempted with that malformed handle.

The initial overbroad expectation that this semantic case must refuse produced
two assertions failures, retained in `probes-discovery.log` and
`probe-summary-discovery.tsv`. The final test distinguishes immutable-bundle
integrity refusal from inherited semantic FAI acceptance and proves baseline
equivalence. Two warnings from a deliberately damaged quoted alias header were
also retained in the discovery log; the final probe suppresses those expected
parser warnings while checking the same refusal.

A separate native hardening change can validate FAI contig uniqueness, positive
line widths, offset/range bounds and FASTA layout consistency with dedicated
tests. This R-only candidate does not change that code or trust boundary. Harvest
additionally pins trusted bundle manifests and artifacts with SHA-256 outside MSS;
the coordinator confirmed public users cannot upload references/start MSS.

## Reproduction and source pins

From the independent evidence directory:

```sh
MSS_TEST_RLIB=/private/tmp/harvest-Rlib Rscript candidate/tests/standalone/preservation.R candidate
MSS_TEST_RLIB=/private/tmp/harvest-Rlib Rscript candidate/tests/standalone/registered.R candidate
Rscript run-probes.R candidate .
```

`run-probes.R` uses local Rcpp/testthat from `/private/tmp/harvest-Rlib`; adjust that
one library location when reproducing elsewhere. The standalone runners accept
`MSS_TEST_RLIB` and otherwise use normal R library discovery. The FAI follow-up
reproduction is contained in the fifth independent test block.

| Source | SHA-256 |
|---|---|
| `R/preservation.R` | `2e2ecf12870189ed0fdd539571eb6ade800ee21f6899a5023c28ac3978407017` |
| `R/projection.R` | `bd367596159973209ef5737edcd1e3e032fb4962fb5869b334207dc5402bc8ca` |
| `src/preservation_core.h` | `9c06526ec9185940d36e3002ff5d9c06a395d34e6a29d1cddb4b7ea119480011` |
| `src/preservation.cpp` | `6d7281661962b2db9f74ae53569841336d810315ab417878b51969d82f63a49e` |
| `src/projection.h` | `766fea8eb745d9f8e18ed49a7091aafd2ac5a935c4249f79e2e53408a8d86997` |
| `src/RcppExports.cpp` | `2b0a3d1bbf1809120dd60d1b1c802b4355f958f2547e545262909feb0db5a1d3` |
| `R/RcppExports.R` | `895594e02ee6c638600abe6085e51d1b0243e706a8dbdccd0d353e9de3b568f2` |

Evidence lives in `results/mss-reference-relocation-independent-review-stage7/`:
`standalone.log`, `registered.log`, `suite-summary.json`, `probes.R`,
`run-probes.R`, `probes-final.log`, `probe-summary.tsv`, discovery logs,
`initial-source-pins.json`, `final-source-integrity.json`,
`author-evidence-verification.json`, `rd-check.log`, and `evidence-index.json`.
No author source, installed library, user repository, active full-reference
preparation or cloud resource was modified. Only tiny review fixtures were copied
or changed. This is a scoped correctness review; no new whole-reference
performance or Windows portability claim is made.
