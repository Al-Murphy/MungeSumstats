# Matched MSS performance benchmark

The timed interval is one complete installed `MungeSumstats::harmonise_file()` call with a verified handle. It includes source hashing and snapshotting, parsing, harmonization, normalized serialization, receipt checksums, and atomic publication. Package loading and verified reference opening have separate elapsed, user, and system times. The Python driver also records process wall time, user and system CPU, and peak RSS across the entire R process. Profiling, external hash comparisons, and oracle work have their own costs and do not enter the processing interval.

The source pins for the first comparison are baseline `a86e8c94d753630eff8c6a695ed2f240d1123d0d` and combined candidate `6abab72a2eb564bb30c6b1d2a98269a7de90a643`. The older deep-profile binary established bottlenecks and planning costs. It does not establish the exact installed baseline for this comparison.

## Before running

Execution requires an admitted, uncontended CPU, memory, and I/O slot. The output volume needs room for all retained editions and scratch files plus the operator’s free-space reserve. Use a new output directory on one filesystem and keep both installations, the dependency libraries, references, and input files immutable for the duration. The driver uses POSIX process high-water counters on Linux and macOS. On Windows select the architecture-specific `Rterm.exe`. The driver uses Win32 process CPU time and peak working set. When a final high-water counter is unavailable it records the observed lower bound and leaves RSS acceptance unverified.

Build both installed package trees from the pinned commits on the same host with the same R, C++17 toolchain, flags, dependency libraries, and install options. Keep those libraries outside the MSS installations. Select the launcher from that R installation for all runs. After build admission, use identical install invocations in the separate source trees and empty package libraries:

```sh
R_LIBS_USER="$COMMON_LIB" "$R_BIN" CMD INSTALL --preclean --library="$BASE_LIB" "$BASE_SOURCE"
R_LIBS_USER="$COMMON_LIB" "$R_BIN" CMD INSTALL --preclean --library="$CAND_LIB" "$CAND_SOURCE"
```

Archive an operator build receipt for each normal, uninstrumented install. Its JSON fields are:

```json
{
  "source_commit": "<full commit>",
  "dll_sha256": "<loaded installed DLL or shared object SHA256>",
  "package_tree_sha256": "<installed package tree SHA256>",
  "toolchain": {
    "compiler": "<compiler executable and version>",
    "compiler_flags": ["<actual CXX17/CXX17FLAGS and relevant CPP/link flags>"],
    "R_version": "<exact R.version.string>",
    "R_platform": "<exact R.version$platform>"
  }
}
```

The driver requires identical `toolchain` objects and checks the observed R version and platform against them. Store sanitized compiler and install transcripts beside the receipts. Record the host CPU, RAM, storage type, local time, and conditions for serial execution in the operator's contention note. Python 3.9+ and the installed R libraries are sufficient for the driver. The driver retains the launcher, actual engine, R home, locale, package versions and namespace locations, numeric libraries, and thread settings.

Prepare the runtime pins with the same launcher used by the worker. On POSIX:

```sh
bench=benchmarks/2026-11-19/performance-upgrade
"$RSCRIPT" --vanilla "$bench/worker.R" runtime "$NEW_RUNTIME_TSV"
```

On Windows use the architecture-specific executable, `Rterm.exe --vanilla --slave --file=worker.R --args runtime NEW_RUNTIME_TSV`. The probe initializes only R and writes a quoted UTF-8 `key,value` TSV. Pin the canonical `R_engine` path as `r_engine`, separately from `r_executable`. Use one `runtime_files` entry with `role,path,sha256` for each available reported path: `r_shared` from `R_shared`, `blas` from `BLAS_library`, and `lapack` from `LAPACK_library`. Omit roles reported as `unavailable`, for example when the numeric library is internal. Linux resolves the actual executable and mapped libR from the R process's own procfs entries. macOS and Windows use the executable reported by R, and associate libR with its installed R home or direct Rterm directory. The probe records that resolution method in `R_shared_resolution`. Every worker checks its engine and library paths against the configured pins before timings enter the report.

From the checkout that contains this directory, after pin and run admission:

```sh
bench=benchmarks/2026-11-19/performance-upgrade
python3 -B "$bench/run.py" pin "$BASE_LIB/MungeSumstats" "$CAND_LIB/MungeSumstats" "$RSCRIPT"
python3 -B "$bench/run.py" pin "$BASE_DLL" "$CAND_DLL" "$BASE_BUILD_RECEIPT" "$CAND_BUILD_RECEIPT"
python3 -B "$bench/run.py" run "$MATRIX_JSON" "$NEW_EVIDENCE_DIR"
```

Set each variable to an explicit local path. The installed tree identity hashes sorted relative file names and their SHA256 digests. `pin` accepts explicit research inputs and each file in the reference inventory after their reads are admitted. The execution driver verifies all declared file hashes and both installed trees before and after all timed processes. It also checks the loaded DLL's path and MD5 from inside each worker. Pre-run pin verification reads the files, and the OS page cache is uncontrolled and never flushed. Fresh processes and handles describe application state. Pair order alternates baseline/candidate, candidate/baseline, baseline/candidate within a three-pair case. Capture operator-observed concurrent work or interference and treat affected timings as such.

Use [matrix.example.json](matrix.example.json) as the settings template. The three scale inputs carry the historical SHA256 identities. Replace path and digest placeholders with admitted local assets, and set the reserve and output estimate for the actual run. Pin `manifest.rds`, `core-version`, `fasta.fai`, `keys`, `offsets`, `records`, `rejected.tsv`, `source.tsv`, the FASTA, and `contig_aliases.tsv` when present. Default opening follows the recorded manifest FASTA and publishes an ordinary receipt. `bind_fasta: true` calls `open_reference(..., fasta=...)` in both revisions and requires the corresponding relocation receipt. Both modes must open the pinned FASTA. Fixture source declarations require the `source_evidence` field and the usual public argument validation. `columns` is the complete supplied mapping when present. Settings default to the public API defaults, with the scale cases explicitly using historical EAF and beta semantics. Projection and remote reference access are outside these timed workloads.

## Workloads and cardinality

| Case | Baseline/candidate pairs | Reference |
| --- | ---: | --- |
| Repeated 1M and 10M chr22 panel | 3 for each size | Format 2 |
| Repeated 40M chr22 panel | 1 | Format 2 |
| Distinct shuffled keys, coordinate-only, ID-only, and negative lookups | At least 1 for each category and format | Formats 2 and 3 |
| Dense contexts, mixed variants, and independent hard cases | At least 1 for each category and format | Formats 2 and 3 |

The repeated input is a bounded panel of 437,458 real chr22 SNP records with synthetic beta and EAF. It reports millions of rows from repeated keys. For the distinct shuffled input use one row per measured eligible key. Its receipt must have equal `distinct_keys` and `input_rows`, and a nonempty `shuffle_derivation` stating the seed or reproducible recipe. For the coordinate-only input omit identifiers. For the ID-only input omit coordinate values while retaining allele, effect, and frequency evidence. For negative cases establish that the tested coordinate and/or ID keys are absent. Preserve the actual distributions of dense sites, mixed alleles, rejection modes, and cache ordering. Format 2 and 3 cases may share one input and independent oracle via `reference_formats` when their expected results agree. Prepare format-specific oracle cases explicitly when warranted.

For each workload provide a pinned input-bound cardinality receipt. The example uses a reference-panel event definition for scale workloads. State independently measured input-key counts and their exact definitions for supplemental cases, including how missing fields, allele orientation, and aliases enter the key. Example:

```json
{
  "input_sha256": "d110828a8b17bf88309db042e7073faa267a0b63a3550d14a87f7761074c91f2",
  "input_rows": 1000000,
  "distinct_keys": 437458,
  "key_definition": "Eligible chr22 panel reference records keyed by chromosome, position, REF, ALT and current ID",
  "method": "<pinned panel inventory and independently verified repeated-input generation receipt>"
}
```

The driver checks input-hash and row-count binding and retains the receipt. The receipt's method supplies the cardinality evidence. Full unique genomewide coverage remains dependent on complete references, representative study inputs, and execution capacity.

Hard-case oracle TSVs have `source_row` first and all 18 other normalized fields in any order, with one ordered record per source row. The driver rejects an incomplete header before starting a timed worker. Write exact expected tokens, including `NA`, signed effects, oriented or unknown EAF, flags, and nullable fields. Cover alias and identifier fanout, repeat indels, multiallelic and ALT–ALT contexts, dense/negative cases, trusted declarations supported by the fixture, specialist preservation, and tricky source/numeric tokens. The oracle must come from hand-derived truth or the independent scientific fixture/verifier, with a retained derivation receipt. Never bootstrap oracle values from either timed installation. Exact full input snapshot and decoded source comparisons also cover raw numeric spelling, quote escaping, CRLF/EOF and gzip where these are present. Run the wider integrity, parsing, API, projection, and native Windows parity suites independently.

For the repeated panel, the tracked independent SNP oracle can check all source and normalized fields from an edition after byte parity. Set `SNP_ORACLE_PROFILE` to a small JSON document with `reference: {source, fasta, fai, build}` and `stages: {1M: {rows, input: {path, sha256}}, 10M: {...}, 40M: {...}}`. Then run one receipt per size under separate verification admission:

```sh
python3 -B inst/verification/fai-safety-stage8/whole-file/harness/verify.py \
  --profile "$SNP_ORACLE_PROFILE" --stage 1M \
  --edition "mss7=$NEW_EVIDENCE_DIR/matched-1M/pair-01-candidate/edition" \
  --receipt "$NEW_EVIDENCE_DIR/snp-oracle-1M.json"
```

## Retention and interpretation

Keep the complete evidence directory. It holds the configuration, pre/post pin inventories, exact worker inputs and argv, raw per-process measurements and logs, R heap diagnostics, all raw editions, stable RDS comparison, full byte digests, audit verification, and every repetition in `summary.json`. `manifest.rds` has a volatile UTC timestamp. All stable RDS fields compare with `identical()`. The returned manifest additionally contains its per-run edition path, which is checked against publication. Edition metadata and reference paths are identical within a case. The driver compares every byte of all other artifacts, independently verifies input and decompressed-source hashes, scans each byte-equivalent normalized edition group for all 19 fields and row ordinals, and checks audit counts against observed statuses.

The 1M/10M target is baseline median processing time divided by candidate median processing time of at least 2, with median candidate peak RSS at or below median baseline RSS. Raw pair ratios and maximum RSS remain in the report. One 40M pair measures scale with a single observation. Hard-case timings need the actual data distribution and timer resolution context. The summary's `target_status` is `pass`, `miss`, or `indeterminate` with exits 0, 2, or 3. A known failed elapsed or RSS criterion establishes a miss even when another criterion is unresolved. Inconclusive RSS or timer resolution otherwise yields an indeterminate result. Successful diagnostic and legacy suites report `not_applicable` and exit 0. Exit 1 indicates a harness or verification failure. Keep valid raw evidence after failures and use a new output directory for any rerun.

Historical completed 1M and 10M editions occupy about 0.234 and 2.345 GB each on a logical-byte basis. The 40M historical case had a 12.228 GB owned-case high-water mark. Six plus six plus two runs give a 39.928 GB (37.19 GiB) scale-only planning figure. Size supplemental panels and logs separately. Budget for scratch, allocated bytes, input/reference copies if needed, independent oracle execution, any diagnostics, and the operator's reserve. The driver checks the provided estimate at startup and the reserve before each R launch. Live disk/RSS caps and contention admission remain the operator's responsibility.

For separate R allocation and CPU diagnostics use a reduced, explicitly documented diagnostic configuration with ordinary paired settings:

```sh
python3 -B "$bench/run.py" run "$DIAGNOSTIC_JSON" "$NEW_ALLOC_DIR" \
  --suite diagnostic --profile alloc --allocation-threshold 65536
python3 -B "$bench/run.py" run "$DIAGNOSTIC_JSON" "$NEW_CPU_DIR" \
  --suite diagnostic --profile cpu
```

`Rprofmem` samples qualifying R allocations and small-object-page markers. `Rprof` samples R execution. Use the standalone `tests/standalone/lookup_parity.cpp` control with candidate `MSS_LOOKUP_COUNTERS` for native search and allocation counters under applicable build and run admission. Counter results measure that fixture's mechanism. If the target is missed, investigate the measured bottleneck using a separately profiled matched workload before changing source again.

Adapted legacy work has its own `--suite legacy` configuration. Supply `legacy_panel` as a pinned biallelic SNP source table with `id, chr, pos, ref, alt` per case and a separate output directory:

```sh
python3 -B "$bench/run.py" run "$LEGACY_JSON" "$NEW_LEGACY_DIR" --suite legacy
```

Both legacy revisions use the historical GRCh37, dbSNP144, one-thread QC settings and the same loader adapter. Panel acquisition remains within the timed `format_sumstats()` call in each revision. The legacy driver checks returned counts and exact emitted table bytes. Legacy output, scientific work, and an isolated chromosome-mapping diagnostic each have separate reporting scopes.
