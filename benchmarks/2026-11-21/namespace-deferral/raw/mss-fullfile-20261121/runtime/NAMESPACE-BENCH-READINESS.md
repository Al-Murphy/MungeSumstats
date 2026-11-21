# Matched Darwin benchmark for VariantAnnotation deferral

Use new normal installs of the ANN-only commit `14ea7af51a19272c7792fbe582a3d42b291b4a1a` and the final committed namespace candidate, both from clean Git archives. Build and measure against `/private/tmp/mss-parallel-support-20261121/common`. The retained inventory records exactly 96 packages. Its 94 inherited package digests match the earlier runner inventory, and the two additions are BiocBaseUtils and GenomicFiles. Keep the existing functional role builds as their own evidence. [Current dependency inventory](/private/tmp/mss-parallel-support-20261121/inventory.json) [Earlier inventory](/private/tmp/mss-fullfile-20261121/fullcall-1m-01/dependency-inventory.json)

## Build and provenance

Let `BENCH=/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade` and `RHOME=/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources`. The frozen `operations/build.py` CLI requires Linux, an image digest and the a86 baseline. Its helpers work on Darwin with real commits. The historical Darwin fresh-install and metadata commands provide the matching template. [Linux CLI constraints](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/operations/build.py:92) [Previous Darwin install command](/private/tmp/mss-fullfile-20261121/runtime/prepare-fullcall.py:60) [Previous Darwin receipt assembly](/private/tmp/mss-fullfile-20261121/runtime/prepare-fullcall.py:135)

For each role create a separate empty library and archive from the same local Git repository only after the candidate commit exists. Record both commits and their tree IDs. Use the same snapshot command and the same normal staged install arguments. Git embeds the archived commit ID, which can be checked with `git get-tar-commit-id < "$archive"`. For direct archive installation, use the package-directory prefix from the historical Darwin builds. The alternative helper `build.archive_source(archive, sha256, commit, new_source_dir)` expects an unprefixed verified Git archive and extracts it for installation from a source directory.

```sh
git -C /private/tmp/mss-lazy-variantannotation-20261121 \
  archive --format=tar --prefix=MungeSumstats/ --output="$archive" "$commit"
/usr/bin/env SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX27.2.sdk \
  "$RHOME/bin/R" CMD INSTALL --use-vanilla --library="$role_lib" "$archive"
```

Use the pinned Command Line Tools and SDK after fresh verification, plus the same clean build environment with `R_LIBS` set to the physical 96-package root, site/user libraries and user startup/Makevars disabled, `MAKEFLAGS=-j1`, and the same locale and thread limits. Check both guard terminals and staged temporary/final package load checks, and compare the actual compiler/link argv and effective flags. Apply the same capture commands once to the common toolchain: compiler version, R `CMD config --all`, and `CC CFLAGS CPPFLAGS CXX17 CXX17FLAGS LDFLAGS`. [Existing compiler check](/private/tmp/mss-performance-20261119/tools/parser-precedence-stage176/build_probe.py:417) [Existing Darwin build receipts](/private/tmp/mss-fullfile-20261121/fullcall-1m-01/candidate-build.json)

Import the frozen `operations/build.py` helpers by putting `$BENCH/operations` on the Python import path. `pin(path)` returns an absolute physical path and SHA-256. `package_tree(role_lib / "MungeSumstats")` returns the runner's strict tree hash and file inventory. `inventory_library(common)` yields the same per-package hash definition as the verified 96-package inventory. It also works on `RHOME/library`. The common library and framework inventories can be revalidated outside the timed phases. For a fresh runtime observation, the existing CLI is `"$RHOME/bin/Rscript" --vanilla "$BENCH/worker.R" runtime "$fresh_tsv"`. The equivalent helper `runtime_roles(new_meta_root, rscript, engine)` runs that command itself and produces the role-labelled shared R, BLAS and LAPACK pins. Cross-check the actual version, platform, R engine and runtime DSO paths against the existing matrix runtime identities. [Build helpers](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/operations/build.py:21) [Runtime identity mode](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R:59)

Write one fresh build receipt per role. The consumed fields are `source_commit`, `dll_sha256`, `package_tree_sha256` and `toolchain`. `toolchain` must be structurally identical for both roles and include `compiler`, `compiler_flags`, `R_version` and `R_platform`. Include actual common compiler/SDK/Makeconf/config pins and dependency inventories there. Keep role-specific tree ID, archive pin, install argv, transcript and guard closure pins alongside the consumed fields. The config entries `baseline` and `candidate` each contain `commit`, `library`, `dll` (pin of `library/MungeSumstats/libs/MungeSumstats.so`), `package_tree_sha256` and `build_receipt` (pin). Keep separate physical role paths. [Exact runner provenance check](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/run.py:232)

## Config subset and launch

Derive a new schema-1 diagnostic JSON from [prep-01/matrix.json](/private/tmp/mss-fullmatrix-20261121/prep-01/matrix.json) by selecting the existing entries below. A `Pin` uses `{"path": "<absolute physical path>", "sha256": "<lowercase hex>"}`.

| Field | Value |
| --- | --- |
| `schema`, `edition_encoding` | Retain `1`, `utf-8` |
| `contention_note` | Fresh admission and observed host/cache conditions for this experiment |
| `r_executable`, `r_engine`, `runtime_files` | Existing runtime pins after fresh matches, or fresh confirmed equivalents |
| `baseline`, `candidate` | The two new role entries above, with their literal Git commit IDs |
| `dependency_libraries`, `dependency_inventory` | Array containing only the physical 96-package common root, and its inventory pin for operator verification |
| `minimum_free_bytes`, `estimated_output_bytes` | Fresh guard-backed reserve and retained-output estimate for this small roster |
| `references` | Retain `matched_2`, `scientific_2`, `scientific_3` verbatim |
| `cases` | Retain exact `hard-cases`, `declared`, and `matched-1M` entries, in that order. Set `pairs: 3` and `target: false` on all three |

With both tiny gold inputs and both reference formats, this expands to five cases and 30 installed workers. Keep `reference_formats`, gold pins, source declaration evidence, cardinality receipts and settings as supplied. The matched file is the fixed 1,000,000-row input `d110828a8b17bf88309db042e7073faa267a0b63a3550d14a87f7761074c91f2` with the format-2 manifest `d15ad358b438fef485ddf0980043556a43f3e4d6b72faca35f3429ce4919ab93` and its original explicit FASTA binding. A prior equivalent six-edition 1M diagnostic retained 1.53 GB peak allocated and used an admitted 3 GiB output allowance and a separate 186 GiB reserve. Use fresh observations for this case. The 42 GB matrix output estimate and warm reserve describe a different roster. [Historical scoped resource receipt](/private/tmp/mss-fullfile-20261121/runtime/MATCHED-1M-RESULT-REVIEW.md:34)

Two explicit bindings are necessary for using the original `run.py` on this experiment. Its CLI hardcodes the a86 baseline even in diagnostic mode. It also requires every worker in a case to have identical loaded namespaces, although the deferred role is expected to have fewer. Launch a small, reviewed input-binding script using `python3 -B <binding.py> run <config.json> <fresh-output> --suite diagnostic --profile none`. It imports the frozen `run` module on its real source path, binds `run.BASELINE` to the reviewed ANN-only commit, and wraps only `run.dependency_state` as below before calling `sys.exit(run.main())`. Record the exact binding source and digest and the resulting role identities. Both config and producer receipts keep their real commit IDs. [Historical baseline assertion](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/run.py:116) [Cross-role namespace assertion](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/run.py:437)

```python
from pathlib import Path
import sys

sys.path.insert(0, "/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade")
import run

original_state = run.dependency_state
eager_only = {
    "AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite",
    "VariantAnnotation", "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "png",
}
observed = {}

def experiment_state(run_dir, config, role):
    case = run_dir.parent.name
    state = original_state(run_dir, config, role)
    framework = Path(config["r_engine"]["path"]).parents[2] / "library"
    library_roots = {framework, *map(Path, config["dependency_libraries"])}
    if any(Path(location) not in {root / name for root in library_roots}
           for name, (_, location) in state.items() if name != "MungeSumstats"):
        raise ValueError(f"Namespace loaded from an unexpected root for {case}/{role}")
    if state != observed.setdefault((case, role), state):
        raise ValueError(f"Namespace identity changed for {case}/{role}")
    baseline = observed[(case, "baseline")]
    if not eager_only <= baseline.keys():
        raise ValueError(f"Expected eager namespaces missing for {case}")
    if role == "candidate":
        omitted = baseline.keys() - state.keys()
        if (omitted != eager_only
                or any(baseline.get(name) != identity for name, identity in state.items())):
            raise ValueError(f"Unexpected namespace delta for {case}: {sorted(omitted)}")
    return baseline

run.BASELINE = "14ea7af51a19272c7792fbe582a3d42b291b4a1a"
run.dependency_state = experiment_state
sys.exit(run.main())
```

The guarded installed namespace probes observed this exact 13-name eager-only difference and zero deferred-only names at both checkpoints. Eager/deferred counts were 87/74 before loading the six references and 93/80 after. Each common namespace had the same version and physical location across roles. MungeSumstats had the same version and resolved from its own verified role library. Each role added exactly the six reference namespaces, with earlier entries unchanged. The probe loaded testthat before the first snapshot and ran against the 102-package union, which exposes the same 96 physical support packages through links and also exposes the six references. Its R code checks each namespace's physical path/version against the union inventory. The frozen public worker uses the 96-package common root directly, loads Rcpp before MungeSumstats, and records its post-call inventory without explicitly loading testthat or the six references. Confirm the exact observed delta in the first tiny worker pair before interpreting public benchmark timings. The three additional observed names beyond the static import-graph difference are `cachem`, `fastmap`, and `memoise`. [Probe source and snapshot boundaries](/private/tmp/mss-reference-closure-20261121/attempt-03/namespaces.R:57) [Probe admission](/private/tmp/mss-reference-closure-20261121/attempt-03/admission.json) [Eager guard receipt](/private/tmp/mss-reference-closure-20261121/attempt-03/eager-namespaces.jsonl) [Deferred guard receipt](/private/tmp/mss-reference-closure-20261121/attempt-03/deferred-namespaces.jsonl)

| Evidence | SHA-256 |
| --- | --- |
| [Namespace probe source](/private/tmp/mss-reference-closure-20261121/attempt-03/namespaces.R) | `279ab8869dc2d35c68a4d9d79fcc8dc7b81ffe087136f5989d2b02b708e89e14` |
| [Eager namespace observations](/private/tmp/mss-reference-closure-20261121/attempt-03/eager-namespaces.json) | `f366b1a7083252c2593baa0da5b525fcc157171e3f21225c4e8e221c8ce80174` |
| [Deferred namespace observations](/private/tmp/mss-reference-closure-20261121/attempt-03/deferred-namespaces.json) | `e452937dc8885bf9175db4ab32af4feb68d074624accdbdcfc3b1785bf8c4fd6` |
| [102-package union inventory](/private/tmp/mss-reference-closure-20261121/attempt-02/union-inventory.json) | `d30550f4c01897764c858eefe03a1d3f037d3cdf0ee90be1224a419c7a2a3213` |
| [Physical 96-package common inventory](/private/tmp/mss-parallel-support-20261121/inventory.json) | `181de043a869c040e467039cd533d936efe3a69734d7c9c27c6a2a80d4e84472` |

The gate checks exact equality across the full version/path inventory within each role and case, with only the role MSS path normalized by the existing parser after identity verification. Across roles it requires the precise eager-only set, rejects any candidate-only namespace and requires a version and physical path match for every shared namespace. Every raw `namespaces.tsv` stays available. If the first tiny pair differs, the gate stops before the later cases and preserves the raw inventories for diagnosis. The frozen worker, its timing boundaries and wait4 accounting are unchanged. [Frozen namespace parser](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/run.py:335) [Public worker setup](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R:99)

## Measurements and interpretation

Run before/after dependency inventory verification and package preparation outside measurements. Finish other test, compile and hashing processes before launch. Run all workers serially on the same hardware with one thread and the existing clean child environment. The frozen driver orders pairs eager/deferred, deferred/eager, eager/deferred in each case. This alternates first position with one extra eager-first pair. Retain pair order and every raw observation along with medians and variability. The OS cache is uncontrolled and the runner pre-reads pinned assets.

Report setup (library paths, both library calls and data.table thread setting), verified reference open, complete public `harmonise_file()` processing, and whole-process wall and Darwin `wait4 ru_maxrss` separately. The tiny processing values can approach timer resolution. Use the per-case entries in `summary.json` and complete peak RSS receipts. Parity requires exit 0, guard closure, normal byte and stable manifest gates across all six editions per case, and exact independent rows for the tiny golds. The matched 1M config has no embedded independent oracle; its frozen checks include all row framing, audit and edition byte parity. Keep any additional pinned SNP oracle outside the timed phase. `target_status=not_applicable` is expected in diagnostic mode. [Frozen worker timing boundaries](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R:99) [Darwin RSS accounting](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/process_metrics.py:101) [Verification and summary](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/run.py:457)

For a separate pure namespace measurement, use fresh guarded `Rscript --vanilla` workers with the same role/common roots and thread settings. Before starting the clock, bind and check `.libPaths()` and capture `loadedNamespaces()` with both MungeSumstats and VariantAnnotation absent. Time only `start <- proc.time(); ns <- loadNamespace("MungeSumstats", lib.loc = role_lib); elapsed <- proc.time() - start`. After the timer, assert the selected namespace and DLL paths and record all namespace names, versions and paths, plus the expected eager/deferred VariantAnnotation state. Pair the invocations in the same alternating order, preserve raw elapsed and process receipts, and keep this phase separate from the unchanged public worker and its broader setup measurement.
