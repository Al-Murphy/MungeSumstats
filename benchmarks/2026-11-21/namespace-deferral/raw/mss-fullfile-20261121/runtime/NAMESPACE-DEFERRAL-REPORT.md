# VariantAnnotation namespace deferral experiment

## Phase 1 regression controls

The experiment is in `/private/tmp/mss-lazy-variantannotation-20261121` on `feat/lazy-variantannotation-20261121` at base `a6de747fb2e7f0bd9a423baa9782c93a6d109795`. Only [the installed regression script](/private/tmp/mss-lazy-variantannotation-20261121/tests/standalone/test_variantannotation_deferral.R) has been added to that worktree.

The `startup` mode checks the fresh MSS load, runs both tabular checks, checks the namespace again, and then uses VCF features without attaching VariantAnnotation. It should fail immediately on the eager BFA install with the exact diagnostic `VariantAnnotation loaded while loading MungeSumstats`. The separate `semantics` process omits that deferral requirement so it can establish a positive reference on BFA. It checks row order, coordinate and effect types, raw source tokens, and the retained invalid row for both reference formats with one-row and single resident batches of frame and data.table input. It also covers custom column mapping and its frequency precedence, real local VCF S4 access and a two-object merge, parsed ANN, indexed input, and an indexed MSS writer round trip with effect and ID metadata. Only a new two-record VCF with ANN is synthesized to drive the affected ANN branch. Repeated tabular operations check values and warnings after VCF loads. The receipt also records the complete normalized and converted frames for exact BFA-to-experiment comparison, along with selected portable provenance and the warning messages.

The current pinned common root lacks GenomicFiles and the two BSgenome reference packages. The installed test script has separate `parallel-vcf` and `bp` modes for the real production parallel reducer and reference-backed BP checker. These are pending those dependencies, as are full legacy formatting with a real reference, SNPlocs and liftover qualification. The ordinary single-threaded VCF path can emit an existing requireNamespace message for GenomicFiles. The S4 split/merge in the core check exercises two real VCF objects through the same registered VariantAnnotation generic. Full integration with MungeSumstats' parallel reducer requires the pending mode. The previous frozen 32-hard/9-declared file golds are separate controls. This script prepares a portable four-row resident fixture to isolate the affected multi-batch bind.

### Proposed baseline commands

Run the two R invocations sequentially under root's existing bounded guard, using a new owner-controlled scratch path for logs and the receipt. The following are the actual R command lines to wrap.

```sh
rscript=/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript
probe=/private/tmp/mss-lazy-variantannotation-20261121/tests/standalone/test_variantannotation_deferral.R
bfa=/private/tmp/mss-fullfile-20261121/builds-01/candidate/library
common=/private/tmp/mss-fullfile-20261121/restored-01/common
scratch=/private/tmp/mss-namespace-phase1-20261121
mkdir -m 700 "$scratch"

env -i PATH=/usr/bin:/bin LC_ALL=C TZ=UTC TMPDIR="$scratch" \
  "$rscript" --vanilla "$probe" startup "$bfa" "$common"

env -i PATH=/usr/bin:/bin LC_ALL=C TZ=UTC TMPDIR="$scratch" \
  "$rscript" --vanilla "$probe" semantics "$bfa" "$common" "$scratch/bfa.rds"
```

The probe resets and verifies normalized installed library roots before requiring MSS. The proposed phase 2 startup command passes both a new candidate receipt path and the BFA receipt path. No R, installs or performance measurements ran during test authoring.

### BFA execution and ANN diagnosis

Root's guarded startup process exited 1 with the expected eager-load diagnostic. The independent semantics process exited 1 while converting the local ANN VCF. Its cbind reported `arguments imply differing number of rows: 9, 2`. The original run and its logs remain the record of this failure. [BFA semantic stderr](/private/tmp/mss-namespace-phase1-20261121/bfa-semantics.stderr)

The source supplies two hypotheses. The inherited VariantAnnotation rbind could assemble the ANN vectors differently from base rbind. Also, the expression `df[, colnames(df) != "ANN"]` uses a computed data.table j expression without `with = FALSE`. If it evaluates as a logical vector, its length is the number of columns, which would explain the 9-versus-2 cbind error. The diagnostic compares both binds with the exact same arguments and inspects the result of both the original and explicit data.table selection. The fixture provides one ANN string with three values for each of its two records and declares three corresponding header fields. The actual imported classes, sizes, labels and values are reported by the diagnostic. It changes no installed code.

```sh
diagnostic=/private/tmp/mss-fullfile-20261121/runtime/namespace-deferral-ann-diagnostic.R
env -i PATH=/usr/bin:/bin LC_ALL=C TZ=UTC TMPDIR="$scratch" \
  "$rscript" --vanilla "$diagnostic" "$bfa" "$common" "$probe"
```

Use the baseline variables above, a fresh root-owned diagnostic scratch directory and the same bounded guard. If the two binds agree on two rows and the original j returns a logical mask, the semantic baseline exposed an existing ANN column-selection bug. Fixing that regression would require a separately approved scope decision with the existing failing ANN check retained. Product source has no edits at this stage.

### ANN prerequisite source change

Root's guarded diagnostic returned raw exit 0 with confirmed closure and a reaped leader (session 59818, terminal 09d1b6). The importer produced two valid ANN rows. Resolved VariantAnnotation and base rbind both returned the same 2-by-3 matrix. The original data.table expression produced a length-nine logical vector. Explicit column selection retained a 2-by-8 data.table and binding the three parsed fields yielded a 2-by-11 data.table. The original binding reproduced the 9-versus-2 failure. [Guarded diagnostic output](/private/tmp/mss-namespace-ann-diagnostic-20261121/diagnostic.stdout)

The prerequisite change at `R/vcf2df.R:118` is:

```diff
-        df <- df[,colnames(df)!="ANN"]
+        df <- df[, colnames(df) != "ANN", with = FALSE]
```

Source SHA-256 for `R/vcf2df.R` at this stage is `4c54077d40071283e4241a8786baa089f2358de73378f7734b4bd1de3d86d295`. The original test script has SHA-256 `761012549d5f3ed81be411752caeb86934078f175b6b4aee547a64ef05ac78b9`. These are the only product and test source paths changed from the experiment base at this stage. The inherited rbind and namespace directives still match BFA. `git diff --check` passed. Root will normally install this prerequisite-only source into a new library and run the unchanged `semantics` mode against it. Its receipt becomes the semantic baseline for the isolated namespace change after the run passes.

Root normally installed the ANN-only build, which passed temporary and final package-load checks (session 44237, terminal 4d4730). The unchanged `semantics` probe passed in its own guarded process with raw exit 0, confirmed closure and a reaped leader (session 49742, terminal 4f9c9e). This produced the first positive semantic baseline that includes the valid ANN fixture. Root authorized local prerequisite commit `14ea7af51a19272c7792fbe582a3d42b291b4a1a` with only the ANN line and the previously authored standalone probe. [Install log](/private/tmp/mss-namespace-annfix-build-20261121/install.stderr) · [Semantic stdout](/private/tmp/mss-namespace-annfix-semantics-20261121/semantics.stdout) · [Semantic guard](/private/tmp/mss-namespace-annfix-semantics-20261121/semantics.jsonl) · [Semantic receipt](/private/tmp/mss-namespace-annfix-semantics-20261121/annfix.rds)

## Phase 2 namespace-only source

The working diff relative to `14ea7af` removes the 11 VariantAnnotation `importFrom` directives and the 10 matching roxygen tags. Ordinary resident batch, header/format mapping, BP and ANN binds explicitly use `base::rbind` (eight call sites). Existing qualified calls, including the production VCF S4 reducer, and the VariantAnnotation installation requirement in DESCRIPTION are intact. Static checks confirmed exact directive/tag symbol correspondence, no remaining bare R rbind calls, identical qualified VariantAnnotation source lines, identical DESCRIPTION and generated help, and a clean `git diff --check` against both base revisions. Roxygen regeneration and installed verification remain pending. [Namespace-only diff](/private/tmp/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-NAMESPACE-ONLY.diff) · [Combined review diff](/private/tmp/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-REVIEW.diff)

The namespace-only diff SHA-256 is `913aa811ef03f50fd1734f3268750f0fb7de99abedcdc25b2b43a4940a44bac8`. The combined review diff SHA-256 is `df9e56442726a2f507d95cb290d45446f4df1c8a89bf46e93c6454be2aab8712`. The namespace changes remain uncommitted pending the installed run and independent review.

For phase 2, root intends matched normal installs of the eager ANN-only prerequisite and the namespace-only candidate against the same newly isolated dependency root including the real GenomicFiles closure. Run prerequisite `semantics` to obtain a matching dependency receipt, then candidate `startup` with that baseline receipt:

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C TZ=UTC TMPDIR="$scratch" \
  "$rscript" --vanilla "$probe" semantics "$ann_fixed_lib" "$common_96" "$scratch/eager-ann.rds"

env -i PATH=/usr/bin:/bin LC_ALL=C TZ=UTC TMPDIR="$scratch" \
  "$rscript" --vanilla "$probe" startup "$candidate_lib" "$common_96" \
  "$scratch/deferred.rds" "$scratch/eager-ann.rds"
```

The same positional output/expected receipt arguments apply to two separate `parallel-vcf` runs and, when the genuine BSgenome packages are available, two `bp` runs. Full scientific/performance and platform evidence is still separate from these probes.

### BP coverage follow-up

The independent source review found no product code issue and requested a BP comparison on both sides of the first VCF call. Root's namespace candidate completed a normal install with raw exit 0 (session 50901, terminal a3d3d2), including its temporary and final package-load checks. The frozen candidate source lives alongside that install. [Candidate install log](/private/tmp/mss-namespace-deferred-build-20261121/install.stderr)

The test-only follow-up records the actual BP results for both GRCh37 and GRCh38 before and after the existing local VCF helper calls the public reader. It asserts identical values, types, row order and warning messages and retains both observations and the intermediate VCF result in the comparison receipt. The reference tests still require the two real BSgenome references. The `semantics`, `startup` and `parallel-vcf` code paths match their first implementation. `git diff --check` passed and the product diff still exactly matches the earlier namespace-only review package. [BP-only diff](/private/tmp/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-BP-COVERAGE.diff)

The probe's new SHA-256 is `81c3eea755faf748b814a3a47672fce511a7901366523885622d97a7d447e9c6`. The BP-only review diff SHA-256 is `4a8f6aecf850bf0fdaf77aa42979e0ca57f776eb3e59a01b73fb2dd27db9e8c6`. Run the updated source probe against the two installed packages with the same real reference closure, first passing `bp "$ann_fixed_lib" "$common_with_references" "$scratch/eager-bp.rds"`, then `bp "$candidate_lib" "$common_with_references" "$scratch/deferred-bp.rds" "$scratch/eager-bp.rds"`. This BP coverage awaits guarded execution.

Root verified the 96-package common dependency inventory, including byte-identical trees for the 94 original packages. Existing installed `test-check_impute_se_beta.R` covers the `es_is_beta = FALSE` format mapping branch. Existing installed `test-read_vcf.R` uses the default `use_params = TRUE`. Root will run both after the real references are ready.

## Darwin runtime acceptance and local candidate commit

Root's guarded runs and independent runtime review passed for the installed eager ANN-fixed and deferred prototypes. The deferred build omitted VariantAnnotation during fresh startup and tabular work and loaded it on the first VCF operation. Core, parallel VCF and reference-backed BP receipts matched byte-for-byte. Each unfiltered installed testthat run passed 269 cases and 3,077 assertions across 89 active files, with zero failures, errors or skips. The two matching existing warnings concerned an SNPlocs import replacement and numeric coercion in the rs-ID test. Extended BP and default VCF/ES mapping coverage passed. The suite used verified copies of the packaged ALS fixture in two affected scripts to keep installed fixtures immutable. Root observed all eight groups and sessions closed. [Independent runtime review](/private/tmp/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-RUNTIME-REVIEW.md) · [Root observation](/private/tmp/mss-namespace-functional-observation-20261121.json)

The authorized local namespace commit is `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b`, tree `f3511124b1fbaf16a7297e92c130e7b521cdf453`, on `feat/lazy-variantannotation-20261121`. Its direct parent is the eager ANN-fixed prerequisite `14ea7af51a19272c7792fbe582a3d42b291b4a1a`. The staged product and BP test changes matched their reviewed patches byte-for-byte. The optional NEWS entry describes first-use dependency loading and the ANN row-preservation fix. The staged set contained exactly those 17 paths and passed `git diff --check`. The branch is clean.

For each role, all 136 committed files under `R/` and `src/` plus DESCRIPTION and NAMESPACE matched the retained prototype source files and preinstall inventory by path, size and SHA-256. Both checkout file sets matched their 1,413-path recorded inventories. The two expected deferred checkout differences recorded in the alignment receipt are the expanded BP probe, which root used from a separately pinned source, and the release note. [File-level alignment receipt](/private/tmp/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-SOURCE-ALIGNMENT.json)

The source-alignment receipt SHA-256 is `fe349f66af494a00a67d9284970b0c420813cb47f44da181dee83bceb7df3494`. The committed candidate is ready for fresh normal builds for matched timing and scientific parity verification. Performance, cross-platform qualification and publication remain pending.
