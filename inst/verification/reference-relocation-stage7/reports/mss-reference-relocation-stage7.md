# Optional verified FASTA relocation — frozen implementation candidate

**Implementation gates GREEN; independent review, full installed-package suite
and native multi-OS CI remain required.** No publication, CI dispatch, adoption,
full-reference copy/build or provider operation was performed by this task.

Frozen commit: **`8334f1c9052ae266dc9df80a35bca07868af1423`**.
Base: `f1cd89cf5e0f24fa754bbe2f0e17c9d93c53d6a6`.
Isolated clean clone: `/private/tmp/gwas-harvest-stage2/mss-reference-relocation-stage7`,
branch `feat/verified-fasta-relocation`. The approved narrow plan is
`mss-reference-relocation-plan-stage7.md`.

## Resulting API and contract

```r
reference <- open_reference("transferred/index", fasta = "transferred/reference.fa")
result <- harmonise_sumstats(input, reference, frequency_type = "EAF")
```

The existing release API had one argument. This appends `fasta=NULL` as the second
argument and preserves existing one-argument behavior; no `verify` parameter is
added or removed. Format 2 and opt-in typed format 3 both work. No basename search,
environment fallback, reference-build inference or partial verification is used.
The recorded original path need not exist when the explicit override is valid.

An explicit override normalizes the physical runtime path, rereads the immutable
manifest even for handle input, verifies the candidate against the original
manifest's exact FASTA size/MD5, verifies every listed FAI/index/source artifact,
then opens fresh native streams. File size/mtime signatures cover the complete
explicit verification/native-open interval. Previous handles, mutated in-memory
metadata and serialized native pointers cannot supply cached verification or a
stale physical FASTA stream to the rebind operation.

The returned handle's `fasta` is its actual normalized physical path. Its optional
`fasta_binding` receipt records:

- `implementation_version="mss-fasta-binding-1"`;
- original `recorded_path` and actual `runtime_path`;
- original verified `size_bytes`, FASTA `md5`, FAI `fai_md5`;
- `verification="full_reference_md5"`.

This is operational provenance. The existing scientific reference checksums,
format, normalizer 7, parser 1, projection 2 and chain MD5 remain unchanged.
The wrapper never rewrites the reference manifest or any bundle payload. The
package uses its existing MD5 integrity contract; Harvest's outer bundle can
independently verify SHA-256 before invoking the package.

Only explicit binding adds `reference_fasta_binding` to table provenance and
streamed manifests, plus `reference_fasta_*` fields in portable edition metadata.
Projection receipts conditionally add `source_fasta_binding` and/or
`target_fasta_binding`, covering both cached handles' actual paths. No new path is
fed into existing scientific reference checksums or chain/projection identities.
Transferring existing editions does not require rewriting their source rows,
normalized rows or historical receipts.

Reusing the opened handle checks its actual file size/mtime even when a subsequent
harmonisation call says `verify=FALSE`. Reopening with an explicit override always
performs the full hash verification again. A worker should therefore open once
and reuse the live handle across batches, as in the existing API.

## TDD and verification actually executed

| Gate | Actual result |
|---|---|
| Initial API RED against untouched release R wrapper | 12 assertion failures, 0 errors; existing default 10 assertions passed |
| First transfer/API GREEN | 7 blocks / 154 assertions |
| Projection-provenance RED before adding optional receipt fields | 8 assertion failures, 0 errors; scientific rows already equal |
| Final relocation blocks inside actual standalone runner | 8 blocks / 184 assertions |
| Final complete standalone preservation/typed/relocation runner | **80 blocks / 1,419 assertions**, no failures/errors/warnings |
| Final generated registered `.Call` runner | **40 blocks / 506 assertions**, no failures/errors/warnings |
| Default old-release R comparison, format 2 and 3 | Complete in-memory results identical; 9 ordinary edition files byte-exact per format |
| Default old-release projection comparison, format 2 and 3 | Complete projection receipts/results identical; 11 projection-edition files byte-exact per format |
| Rd validation and `git diff --check` | Passed |

The reference bundle plus tiny FASTA is moved to a path containing spaces and
UTF-8, with the old location absent. No symlink or container permission is needed.
The tests compare exact source tokens, all 19 native normalized fields, IDs,
effects, nullable multiallelic EAF, unchanged checksums/algorithm identities and
unchanged hashes of **every bundle file including the RDS manifest**. Independent
literal effects are also asserted for the SNP/swapped/indel fixture.

Negative cases include a wrong same-length FASTA, missing file, empty/NA/vector/
noncharacter path and directory, malformed FAI, changed records/keys/offsets/
source, and unsupported reference version. Rebinding a poisoned cached handle
rereads the real manifest; rebinding a serialized handle creates a fresh native
pointer. A changed runtime file signature invalidates the handle even with
`verify=FALSE`.

The projection fixture transfers both GRCh37 and GRCh38 tiny references, tests
both mixed 2→3 and 3→2 format directions, and compares the **entire** projected
normalized table and source. It checks source/target actual handle paths,
unchanged old receipt fields, both new optional binding receipts, the published
projection receipt, and byte-exact normalized/source artifacts.

Default compatibility uses the unchanged release R functions in a separate
environment with the same unchanged native core. Ordinary edition comparisons
cover `normalized.tsv`, `source.tsv`, `source-artifact`, `audit.tsv`, `metadata.tsv`,
`columns.tsv`, `reference_checksums.tsv`, `reference_contigs.tsv`, `checksums.tsv`.
Projection comparison additionally covers `projection.tsv` and
`projection_audit.tsv`. Manifest names/values are identical after removing only
the operational creation timestamp and returned output location. No blanket
claim of RDS byte identity across different creation times is made.

`tests/standalone/preservation.R`, which the actual portable GitHub workflow
executes, explicitly loads the relocation tests after typed-format tests. The
registered runner also loads both suites. No workflow was dispatched here.
The exact registration runner is `tests/standalone/registered.R`: it copies
`preservation.cpp`, `preservation_core.h`, `projection.h` and `RcppExports.cpp`
to a fresh temporary directory, invokes `R CMD SHLIB` with Rcpp headers/C++17,
loads the resulting native DLL, resolves all seven registered symbols, and uses
the actual `R/RcppExports.R`. It requires only Rcpp/testthat and the ordinary R
compiler toolchain (Rtools on Windows), with no Bioconductor load, bcftools, WSL
or container. The coordinator can add a separate portable workflow step using
`shell: Rscript {0}` and `run: source("tests/standalone/registered.R")` after the
existing dependency installation. That is a proposed separate CI coverage change,
not an already executed Windows proof or an edit to this frozen candidate.
Full installed-package validation was deliberately left to the coordinator's
post-review gate, as requested. Windows behavior is represented in portable
fixtures but has not yet been established by an actual Windows run for this
candidate.

## Source identities and scope

| Changed artifact | SHA-256 |
|---|---|
| `R/preservation.R` | `2e2ecf12870189ed0fdd539571eb6ade800ee21f6899a5023c28ac3978407017` |
| `R/projection.R` | `bd367596159973209ef5737edcd1e3e032fb4962fb5869b334207dc5402bc8ca` |
| `tests/testthat/test-reference-relocation.R` | `eda1be518d38efecec0e54c73803164ae75ac7a8bef943d446968675a4bedac5` |
| `tests/standalone/preservation.R` | `1e5b2733c6307b6ae6c27fbb8d70914916e95b111e41e18272b5fee36508dce6` |
| `tests/standalone/registered.R` | `31b36a26d8747881ad329ed894138cef6bbaf7c28f4e78e66e82c47305e7d33a` |
| `inst/PRESERVATION_CONTRACT.md` | `1dbd69cdd9746119acfbc6e699752e248b000941db100db6df41a38a4b66ef1c` |
| `man/open_reference.Rd` | `a73ce0488b370671121111bb09a558d29ecf6f9377b05b42e41a19b62809da95` |

All C++ source and generated native/R registration wrappers are byte-identical
to the base. In particular the core stays
`9c06526ec9185940d36e3002ff5d9c06a395d34e6a29d1cddb4b7ea119480011`;
native projection stays
`766fea8eb745d9f8e18ed49a7091aafd2ac5a935c4249f79e2e53408a8d86997`.
No reference writer, marker, default format, normalization/sign rule or variant
retention behavior changes. Frozen installed libraries and active full reference
preparations were untouched. Only tiny reference fixtures were moved or copied.

Evidence root: `results/mss-reference-relocation-stage7/`. Final files include
`source-pins-final.json`, `test-summary-final.json`, `candidate.patch`,
`standalone-final.log`, `registered-final.log`, both RED logs,
`default-compatibility-final.log/.tsv`, `Rd-check-final.log` and
`evidence-index.json`. Historical working pins/logs precede the final projection
addition and are superseded by `source-pins-final.json` and final successful logs.

The frozen candidate is ready for independent review and the coordinator's
installed-package/native-CI gates. Harvest worker/config integration remains a
separate coordinator task; this implementation changes MSS only.
