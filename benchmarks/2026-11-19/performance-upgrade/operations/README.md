# Reproducing the matched Linux benchmark

Use a dedicated, admitted Linux amd64 CPU host. The November 19 capacity plan quotes Runpod Secure Cloud `cpu3m` in `US-MD-1`, four vCPUs, 32 GB RAM and 100 GB scratch at approximately $4.21 for 18 hours. Use an immutable image digest. The catalog's `runpod/base:1.0.2-ubuntu2404` is a starting image and lacks the prepared R environment. Provisioning, scientific transfer and timed work each require the corresponding resource slot. These helpers perform no provider operation.

## Freeze the inputs and execution environment

1. On an image builder, start from an approved Linux/amd64 base digest with Python 3.12+, R and matching Bioconductor system prerequisites, C++17, GNU make, git, rsync, coreutils and all system development libraries for the package DESCRIPTION. Record the image digest, system package versions and actual compiler. The image may already contain the shared dependency library. Otherwise freeze and install it once before the timed host is allocated. Include the packages required by the full package tests, including their reference data. The shared R library must exclude MungeSumstats. Source biological data-package archives, system libraries and install logs can be much larger than the benchmark inputs. Size them on the builder and keep the archive cache off the 100 GB scratch disk.

   ```bash
   # OPS and BASE_SOURCE are absolute paths on the separately admitted builder.
   # REPOS is a reviewed text file with the chosen versioned CRAN snapshot and
   # matching Bioconductor software, annotation, experiment and workflow repo
   # URLs, one public HTTPS repo per line. Freeze the R/Bioc release in the image.
   Rscript --vanilla "$OPS/dependencies.R" freeze "$BASE_SOURCE/DESCRIPTION" "$DEP_ARCHIVES" "$REPOS"
   sha256sum "$DEP_ARCHIVES/lock.tsv" > "$BUILDER_EVIDENCE/deps-lock.sha256"
   Rscript --vanilla "$OPS/dependencies.R" install "$BASE_SOURCE/DESCRIPTION" "$DEP_ARCHIVES" /opt/mss-deps
   ```

   Retain the repos list, lock, archives, complete per-package logs, build recipe, image digest and clean image-builder logs. `dependencies.R` resolves transitive hard dependencies of Imports/LinkingTo/Depends and all direct Suggests. It replays in topological order after checking every source SHA. Base/recommended packages actually present in the base R library come from the immutable image. Verify the realized image contains the full test data or explicitly account for any missing external fixture before performance admission. Supply the image digest observed by the operator as `IMAGE`; an environment variable alone cannot attest the scheduler actually launched it.

2. Create two uncompressed `git archive` tar files from full Git object IDs. The baseline is `a86e8c94d753630eff8c6a695ed2f240d1123d0d`. The current combined production/test candidate is `f3a5c5cf6ff38920ed07b0cf4058dae61f9dd06e`. Pass the final approved candidate explicitly when generating archives and launching the workload. The benchmark harness itself can have a separate commit. Transfer that exact reviewed checkout, or its committed archive, containing the runner, these helpers and the tracked stage-8 profile/oracle. Pin its full commit and archive digest in the external launch record. Use a directory containing only source, without .git or credentials.

   ```bash
   # Staging system, after archive/transfer admission. REPO is the source checkout.
   set -eu
   BASE=a86e8c94d753630eff8c6a695ed2f240d1123d0d
   CANDIDATE=f3a5c5cf6ff38920ed07b0cf4058dae61f9dd06e
   for rev in "$BASE" "$CANDIDATE"; do
       git -C "$REPO" archive --format=tar --output="$ARCHIVES/$rev.tar" "$rev"
   done
   (cd "$ARCHIVES" && shasum -a 256 "$BASE.tar" "$CANDIDATE.tar" > SHA256SUMS)
   git -C "$HARNESS_REPO" archive --format=tar --output="$STAGE/harness.tar" "$HARNESS_COMMIT"
   shasum -a 256 "$STAGE/harness.tar" > "$STAGE/harness.sha256"
   ```

3. Transfer just the following inventory from the completed stage-2 archive. On the Linux scratch volume use the paths in the right column, with the same structure on export. Use the historical SHA256 values in `inst/verification/fai-safety-stage8/whole-file/profile.json`. Keep source checksums independently at the receiving end.

   | Staging source | Linux path relative to scratch |
   | --- | --- |
   | `/private/tmp/gwas-harvest-stage2/assets/mss/matched-{1M,10M,40M}.tsv` | `inputs/mss/` with the three original names |
   | `/private/tmp/gwas-harvest-stage2/assets/reference/GRCh37.ensembl104.chr22.fa` | `inputs/reference/GRCh37.ensembl104.chr22.fa` |
   | `/private/tmp/gwas-harvest-stage2/assets/reference/matched-index/{manifest.rds,core-version,fasta.fai,keys,offsets,records,rejected.tsv,source.tsv}` | `inputs/reference/matched-index/` with the eight original names |
   | The two package tar files and `SHA256SUMS` | `archives/` |
   | Reviewed harness tar and checksum record | Dedicated harness checkout outside results |
   | Written UTC admission, actual CPU, RAM and filesystem/available bytes, workload/I/O isolation and timing/transfer schedule | `admission.txt` |

   The historical FASTA is 52,159,647 bytes. Copy it byte for byte. The transferred format-2 manifest has the original absolute FASTA path. The operations helper verifies the historical pins, calls `open_reference(directory, fasta=relocated_path)`, and sets `bind_fasta=true` for **both** timed revisions. Preserve the manifest, FAI and every byte of the index. The fresh baseline builds matched format 3 from the pinned format-2 source, FASTA, build and version using its public `prepare_reference`. Its manifest records its new real path, so it uses the standard default open. Both scientific indexes use the same public baseline function.

   Example staging commands from the admitted transfer machine. `HOST` is the already verified SSH host alias of this allocation. `STAGE` is a new local directory; `ARCHIVES=$STAGE/archives`. On the destination provision `/workspace/mss-task4` and `/workspace/mss-harness` as new directories. Keep the reviewed manifest checksum and source commits in the independent launch record.

   ```bash
   assets=/private/tmp/gwas-harvest-stage2/assets
   mkdir -p "$STAGE/inputs/mss" "$STAGE/inputs/reference/matched-index"
   cp "$assets/mss/"matched-{1M,10M,40M}.tsv "$STAGE/inputs/mss/"
   cp "$assets/reference/GRCh37.ensembl104.chr22.fa" "$STAGE/inputs/reference/"
   cp "$assets/reference/matched-index/"{manifest.rds,core-version,fasta.fai,keys,offsets,records,rejected.tsv,source.tsv} "$STAGE/inputs/reference/matched-index/"
   (cd "$STAGE" && shasum -a 256 harness.tar > harness.sha256)
   rsync -a --partial "$STAGE/inputs" "$STAGE/archives" "$STAGE/admission.txt" "$HOST:/workspace/mss-task4/"
   rsync -a --partial "$STAGE/harness.tar" "$STAGE/harness.sha256" "$HOST:/workspace/mss-harness/"
   # On the destination:
   cd /workspace/mss-harness
   sha256sum -c harness.sha256
   test "$(git get-tar-commit-id < harness.tar)" = "$HARNESS_COMMIT"
   tar -xf harness.tar
   HARNESS=/workspace/mss-harness
   ```

## Execute and preserve all phases

Start an external 18-hour deadline at allocation, with a 16-hour workload cutoff and the final two hours for verified export and teardown. Track the exact newly allocated Pod ID, storage and end time in the supervisor's record. Make a clean dedicated scratch directory, a private empty home and private temp directory on that volume. Use direct SSH/rsync or another separately approved transfer method. Verify the archive checksums, historical inputs and checkout identity after transfer. The clean home has no SSH keys, PATs, repository auth or other user configuration. A credential-free workload needs no token environment variables. Keep shell tracing off.

```bash
# Linux. s must be an absolute dedicated scratch mount/directory.
set -eu
umask 077
s=/workspace/mss-task4
ops="$HARNESS/benchmarks/2026-11-19/performance-upgrade/operations"
mkdir "$s/empty-home" "$s/tmp"
(cd "$s/archives" && sha256sum -c SHA256SUMS)
# ALLOCATION_UTC is the actual allocation start with zone, not login time.
# IMAGE is the separately verified sha256:... digest of the provisioned image.
env -i PATH="$PATH" HOME="$s/empty-home" TMPDIR="$s/tmp" LANG=C.UTF-8 LC_ALL=C.UTF-8 TZ=UTC \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1 \
    python3 -B "$ops/guard.py" --scratch "$s" --allocation-utc "$ALLOCATION_UTC" \
    --log "$s/guard.jsonl" -- bash "$ops/sequence.sh" "$s" /opt/mss-deps "$CANDIDATE" "$IMAGE"
```

`sequence.sh` includes all argv for every phase. It requires a new results directory and returns the matrix outcome after the scientific checks and diagnostics succeed. Both package installs run `R CMD INSTALL --preclean` into separate new libraries. Build evidence includes the git archive ID and SHA256, exact R and compiler path/version/flags, successful install argv and full transcript, binary and installed package-tree identities, shared dependency identities, system package and R package versions. Failed install transcripts stay in place. The runner takes final R engine, R shared library, BLAS and LAPACK pins from the revised `worker.R runtime` interface. Use the accepted matching revision of the runner for launch.

| Order | Contents of results | Gate |
| --- | --- | --- |
| 1 | `build/`, `small/`, `smoke/` | Tiny unprofiled pairs for generic/declared fixtures in both formats, known all-19-field gold, pinned runtime and allocated reference, full artifact comparisons |
| 2 | `functional/` | Both revisions' complete installed testthat suites, standalone/registered APIs, native suites, error parity for corruption, parser limits, bad gzip, acquisition, fanout caps and malformed fields, and exact differential lookup output. Each installed suite writes per-test counts, condition messages and its qualification decision. A skip, error, failed expectation or suite without passing assertions stops the sequence before full asset preparation and timing. Warnings remain in the per-role evidence. |
| 3 | `assets/`, `matrix/` | Seven matched scale pairs; four chr22 supplemental inputs and four scientific inputs each in both formats. All paired editions and comparison receipts retained; cross-format normalized bytes checked separately. |
| 4 | `snp-oracle-{1M,10M,40M}.json` | Tracked independent SNP verifier checks the full pinned input, complete FASTA-backed panel and all 19 fields for every row of one candidate edition per size. The runner checks byte identity across all other repetitions. |
| 5 | `legacy/`, `alloc/`, `cpu/`, `lookup-pressure.log` | One separate paired adapted `format_sumstats` 1M comparison, two separately profiled paired 1M diagnostics, and native lookup pressure counters. Timed preservation targets use their original full matrix repetitions. |

The projected retained scale bytes are 39,927,958,508 before supplemental outputs. Allow 49 GB for the whole unprofiled matrix, about 2.55 GB for matched input/FASTA/index transfer, under 1 GB for generated assets, roughly 3–6 GB for diagnostics and retained editions outside the matrix, plus filesystem allocation overhead, build logs, test scratch and at least 8 GB reserve. Source archives and the pinned dependency image add their actual footprint. The previous macOS baseline scale work suggests about 26 minutes for the baseline scale processing intervals alone. For planning allow 1–3 hours for matched builds/qualification on a prepared image, 2–7 hours for matrix processing and verification, 1–4 hours for full SNP oracles, and up to 2 hours for transfer/setup/diagnostics. These are allowances for serial work under the 16-hour cutoff, unmeasured on this host. Heavy dependency provisioning belongs on the separately sized image builder. Use observed throughput/disk/RSS to assess whether capacity fits the complete job.

`guard.py` samples free bytes and owned process-group RSS every second, issues TERM then KILL on failure, and stops at the absolute compute cutoff. It sets a 26 GB address-space ceiling per process, an observed 24 GB group RSS cap and an 8 GB free-space reserve. Sampling permits transient overshoot; the dedicated host's storage and memory limits provide the outer bounds. Use only foreground commands in this group. An out-of-band watcher must enforce the 18-hour limit if the host or in-container guard disappears. Once timed processes have stopped or completed, export the data, input/source pins, complete worker editions, all repetitions, manifests, logs and guards to the approved destination. Between timed phases incremental export is possible; avoid competing transfers during measured intervals. Generate an inventory of the final retained paths/SHA256 at source, then verify at destination. Keep both copies until verification or the hard deadline. The external supervisor tears down only the recorded owned allocation under its authorization and confirms provider state/billing. Container exit does not stop billing. If interrupted, report the exact completed phases, keep all partial evidence and use a fresh results path on any rerun.

```bash
# On the workload host after all writers have quiesced. Move any intentional
# temp/home evidence into results before sealing; keep the manifest separate.
set -euo pipefail
cd /workspace/mss-task4
python3 -B "$ops/retained.py" seal . RETAINED.json
```

```bash
# Receiving machine, identical pinned operations helper, fresh approved export.
set -euo pipefail
mkdir "$EXPORT/task4"
python3 -B "$ops/retained.py" roots > "$EXPORT/retained-roots.nul"
rsync -a --partial "$HOST:/workspace/mss-task4/RETAINED.json" "$EXPORT/RETAINED.json"
rsync -a --recursive --partial --from0 --files-from="$EXPORT/retained-roots.nul" \
    "$HOST:/workspace/mss-task4/" "$EXPORT/task4/"
rsync -a --partial "$HOST:/workspace/mss-harness/harness.tar" "$HOST:/workspace/mss-harness/harness.sha256" "$EXPORT/"
python3 -B "$ops/retained.py" verify "$EXPORT/task4" "$EXPORT/RETAINED.json" > "$EXPORT/task4-verification.log"
(cd "$EXPORT" && shasum -a 256 -c harness.sha256 > "$EXPORT/harness-verification.log")
printf '%s\n' 'Verified retained paths, file contents and harness archive.' > "$EXPORT/VERIFIED"
```

The same five roots define the source inventory and the rsync selection. The receiver independently compares every relative path and object type, including empty directories, before checking file digests. Extra top-level or nested paths, missing paths, symlinks, special files, unreadable entries and changed content fail verification. `RETAINED.json` is written only after reading the entire retained source tree. An interrupted transfer can be retried against the sealed, quiescent source. Keep failed verification logs and the source evidence.

## Fixture derivations

The scale inputs keep the historical byte hashes in the matrix template and stage-8 profile. For each size the generator independently reconstructs the complete source row stream from all 437,458 reference events in reference order, alternating direct and swapped rows at each cycle. A partial final cycle is permitted. The receipt binds exact bytes, source rows and unique canonical events independent of input orientation. All supplemental row counts and pairs stay as specified in the template.

The four supplemental chr22 workloads use one ordered set of eligible source events. Ascending SHA256 over UTF-8 seed `matched-chr22-20261119`, tab and current ID gives a reproducible shuffled order, with ID breaking digest ties. Coordinate-only drops the SNP column, ID-only drops both coordinate columns. The negative panel uses increasing coordinates absent from the complete source-locus set and unique synthetic identifiers checked against all current IDs and aliases. Input-bound receipts count the actual distinct source keys after parsing each file.

`scientific.py` is the independent specification of the synthetic reference and source/gold pairs. At S:2 and S:6 the same alias names two distinct C>T events. A coordinate distinguishes them. At T:1 a single reference A has alternatives C and G, so swapping suppresses EAF and the C/G contrast is ALT–ALT. Q:1 adds T, introducing a competing reverse-strand ALT–ALT explanation. At I:4 and I:5 equivalent AA>A deletions in the A run left-align to I:1 CA>C. The declared fixture actually supplies REF, ALT and strand, including reverse and unknown cases. These examples follow the reasoning in `test-preservation.R`, `test-alt-alt-orientation.R`, `test-source-reference-contract.R` and `test-equivalent-indel-context.R`. All generic and declared expectations are literal, including candidate counts, effects, EAF and flags. Gold construction only expands nulls, source ordinals, entity identity, constant build and four null unprojected fields.

The generic hard input includes 32 literal rows with specialist, numeric overflow/underflow, malformed position, quoted UTF-8 IDs, quoted tabs, CRLF, gzip and missing final LF. The nine extra declared rows are separately tagged and checked in both formats. Mixed cycles the 32 generic rows 3,125 times. Dense cycles four allele orders and 98 identifiers at a single site whose 96 duplicate/split source records yield three canonical ALT events. The observed census includes original contig/ID aliases and ordered alleles. Neither scientific gold nor repeated-panel gold uses outputs of a timed binary.

## Windows qualification

Run native Windows x64 with the same two full source IDs in separate fresh R libraries, matched 64-bit R/Rtools C++17 and common pinned Windows dependencies. Capture the OS/R/compiler/toolchain, source and DLL identities with that host's paths and logs. Use the repository's native `tests/standalone/preservation.R` and `registered.R` commands and the installed package test suite for each revision. Run the identical portable lookup parity fixture against both headers under Rtools and byte-compare the entire outputs. Include the UTF-8/alternate encoding, newline, FAI, relocation, quoted-stream, error/limit and projection cases in the Windows review. For portable runner validation materialize the small fixture matrix in native paths with the correct Windows launcher/engine and runtime roles, then run the diagnostic smoke using `Rterm.exe`. Retain observed peak working set or its explicit incomplete status and file encoding. Linux timings and Windows portability receipts are distinct deliverables. Full genomewide claims require separate data and validation.

```powershell
# Native PowerShell, R/Rtools on PATH; R and Rtools paths/pins recorded externally.
# $sourceRoot has baseline and candidate source subdirectories. $e is new.
$env:R_LIBS_SITE = $deps
foreach ($role in 'baseline','candidate') {
    $src = Join-Path $sourceRoot $role
    $lib = New-Item -ItemType Directory -Path (Join-Path $e "$role-library")
    $env:R_LIBS_USER = "$($lib.FullName);$deps"
    & R.exe CMD INSTALL --preclean "--library=$($lib.FullName)" $src *> "$e/$role-install.log"
    if ($LASTEXITCODE) { throw "$role install" }
    foreach ($suite in 'preservation','registered') {
        & Rscript.exe --vanilla "$src/tests/standalone/$suite.R" $src *> "$e/$role-$suite.log"
        if ($LASTEXITCODE) { throw "$role $suite" }
    }
    & Rscript.exe --vanilla -e 'a<-commandArgs(TRUE);testthat::test_dir(file.path(a[1],"tests/testthat"),package="MungeSumstats",load_package="installed",stop_on_failure=TRUE)' $src *> "$e/$role-testthat.log"
    if ($LASTEXITCODE) { throw "$role testthat" }
    $flags = @('-std=c++17','-O1','-I',"$src/src")
    if ($role -eq 'candidate') { $flags += '-DMSS_LOOKUP_COUNTERS' }
    & g++.exe @flags "$sourceRoot/candidate/tests/standalone/lookup_parity.cpp" -o "$e/$role-lookup.exe" *> "$e/$role-compile.log"
    if ($LASTEXITCODE) { throw "$role lookup build" }
    $p = Start-Process "$e/$role-lookup.exe" -ArgumentList 'parity',("`"" + "$e/$role-fixture" + "`"") -Wait -PassThru -NoNewWindow -RedirectStandardOutput "$e/$role-parity.out" -RedirectStandardError "$e/$role-parity.err"
    if ($p.ExitCode) { throw "$role parity" }
}
if ((Get-FileHash "$e/baseline-parity.out").Hash -ne (Get-FileHash "$e/candidate-parity.out").Hash) { throw "Windows parity" }
py -3 -B "$bench/run.py" run $windowsSmokeJson "$e/runner-smoke" --suite diagnostic
if ($LASTEXITCODE) { throw "Windows runner smoke" }
```
