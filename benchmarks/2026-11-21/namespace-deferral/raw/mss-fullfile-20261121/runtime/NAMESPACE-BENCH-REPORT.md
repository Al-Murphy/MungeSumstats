# Darwin namespace benchmark binding

The binding prepares and launches a diagnostic for the actual ANN-only baseline `14ea7af51a19272c7792fbe582a3d42b291b4a1a` and a distinct committed namespace candidate. Its five expanded cases are `hard-cases-f2`, `hard-cases-f3`, `declared-f2`, `declared-f3`, and `matched-1M`. Each case has three balanced pairs and six fresh workers. The 30 workers retain the matrix's literal tiny golds, input and reference pins, source declarations, cardinalities and call settings. The format-2 matched case retains its explicit FASTA binding. [Readiness audit](NAMESPACE-BENCH-READINESS.md)

The binding authenticates `warm-current.py` at SHA-256 `c58b14cda2ad650cefa078e4a73cd2a2a4d5ab8ffc4f7477e5375e0aaa794220` before importing its verified bytes, then calls its eleven-source authentication and loader. The loaded Python modules retain their physical `__file__` paths. During the diagnostic, the two frozen runner bindings are the ANN baseline ID and the namespace comparison callback. The original parser still verifies the selected MungeSumstats role. The callback also verifies each raw version and physical path against the installed DESCRIPTION of the matching package, retains exact per-case and role inventories across repetitions, and requires exactly the 13 eager-only packages and identical shared metadata. The first deviating worker fails the gate. The comparison view contains only shared namespaces and normalizes the already verified role-specific MungeSumstats path. [Binding](namespace-bench.py)

## Input contract

`PROVENANCE.json` is a root-produced control from actual completed build and runtime phases. A pin is `{"path": "<absolute physical file>", "sha256": "<actual lowercase SHA-256>"}`. Supply these top-level fields.

| Field | Required value |
| --- | --- |
| `repository` | Physical local Git repository containing both actual commits and trees |
| `baseline`, `candidate` | Pins of the separate role build receipt JSON files |
| `runtime_receipt` | Pin of the fresh frozen-worker runtime TSV |
| `runtime_transcripts` | Pins of the actual runtime guard JSONL, stdout and stderr, in that order |
| `contention_note` | Fresh host admission, concurrency and cache observations for this experiment |
| `minimum_free_bytes`, `estimated_output_bytes` | Actual operator-admitted reserve and retained output estimate in bytes. The plan is 186 GiB reserve and 3 GiB retained output. |

Each role receipt has the fields below. [Historical Darwin receipt shape](/private/tmp/mss-fullfile-20261121/fullcall-1m-01/candidate-build.json)

| Field | Required value |
| --- | --- |
| `source_commit`, `source_tree` | Actual 40-hex Git commit and Git tree IDs. The binding checks both in the local repository and checks the archived commit ID. |
| `archive`, `library` | Pin of the actual guarded install's Git archive, and physical path of the installed library containing only `MungeSumstats` |
| `dll`, `dll_sha256` | Pin of `library/MungeSumstats/libs/MungeSumstats.so`, and its SHA-256 |
| `package_tree_sha256` | Actual `operations/build.py` `package_tree(library / "MungeSumstats")["tree_sha256"]` |
| `install_argv`, `install_transcripts`, `terminal` | Exact guarded install argv; pins of guard JSONL, stdout and stderr in that order; actual guard terminal object |
| `compile_and_link_argv` | Three tokenized `clang++` commands from the actual install stdout, matching across roles |
| `toolchain` | One structurally identical object across role receipts with the fields below |

Shared `toolchain` fields are `R_version`, `R_platform`, `compiler` (the compiler version stdout as a list of lines), `compiler_flags` (name-to-value mapping, including `CC`, `CFLAGS`, `CPPFLAGS`, `CXX17`, `CXX17FLAGS`, `LDFLAGS`), `sdk` (physical path), `build_toolchain_pins` (physical common toolchain files including the existing pinned compiler, SDK and R Makeconf), `compiler_receipt` (pin of actual compiler version stdout) and `r_config_receipt` (pin of actual `R CMD config --all` stdout). The binding verifies those flags against the actual capture and compares the runtime TSV to the pinned R, BLAS and LAPACK identities. The pinned SDK must be `/Library/Developer/CommandLineTools/SDKs/MacOSX27.2.sdk` to match the existing toolchain file pin. The frozen worker's native measurement is Darwin `wait4 ru_maxrss`. [Historical toolchain pins](/private/tmp/mss-fullfile-20261121/fullcall-1m-01/dependency-inventory.json)

Use guard labels `normal-install-baseline` and `normal-install-candidate`. Each normal staged install uses `/usr/bin/env SDKROOT=… <RHOME>/bin/R CMD INSTALL --use-vanilla --library=<physical-role-library> <physical-git-archive>` with the same optional historical `/usr/bin/time -l` prefix. Use the same clean build environment with the physical 96-package common directory in `R_LIBS`. The guard records environment keys, and the binding verifies matching keys and an actual Rcpp include path in both compiles. The root operator provides the recorded environment values. The binding checks both committed archive identities, actual role DLL and package-tree digests, equal native compile/link argv, both staged load markers and successful guard closure. Use guard label `runtime` for the direct `Rscript --vanilla <frozen worker.R> runtime <fresh TSV>` invocation.

## Preparation and measurement commands

After obtaining the actual evidence, use the physical pinned Python. `META` is new or empty, private and physical. `RUN_PARENT` is an existing, separate private physical directory. The measured output itself must be absent. Root closes other heavy phases, admits resources, and launches using the existing guard before measurement.

```sh
PY=/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15
BINDING=/private/tmp/mss-fullfile-20261121/runtime/namespace-bench.py
"$PY" -I -B -S "$BINDING" prepare "$PROVENANCE" "$META" --candidate "$CANDIDATE_COMMIT"
"$PY" -I -B -S "$BINDING" run "$META/diagnostic.json" "$RUN_PARENT/diagnostic" --suite diagnostic --profile none
```

Preparation writes `diagnostic.json`, `binding.json` and `preparation-inventory.json`. The binding seals the exact config, provenance and preparation inventory, checks all of them again at execution, and takes before/after inventories of both role installations, the 96-package common library and the 30-package R framework library. The two dependency inventory controls preserve the 94 inherited package digests and the two real additions. [Common inventory](/private/tmp/mss-parallel-support-20261121/inventory.json)

After runner success and the parity, RSS and inventory checks, `namespace-policy.json` records each case's actual loaded role states, the measured delta, the common comparison state and all 30 raw `namespaces.tsv` pins in execution order. Worker TSVs stay in the runner directories. `namespace-inputs-before.json` and `namespace-inputs-after.json` retain the inventory comparison. The guard and operator still close the outer execution separately. A diagnostic success is `target_status=not_applicable`. The independent matched 1M oracle belongs after timing. The full C9 matrix, warm and platform gates stay open.

## Source verification

| Artifact | SHA-256 |
| --- | --- |
| [Binding](namespace-bench.py) | `88b9f264c26c4e0d6c4f84a09ebf9384132c328c36358b66548bdbe6fb0cd998` |
| [Source controls](namespace-bench-source-controls.py) | `a55dad3edbad7b680387d7178d18fb927f52ba5ecaa664b954e3513b48288733` |

```sh
RUNTIME=/private/tmp/mss-fullfile-20261121/runtime
"$PY" -I -B -S "$RUNTIME/namespace-bench-source-controls.py"
"$PY" -I -B -S -c 'import ast, sys; from pathlib import Path; [ast.parse(Path(p).read_bytes(), filename=p) for p in sys.argv[1:]]' "$BINDING" "$RUNTIME/namespace-bench-source-controls.py"
"$PY" -I -B -S "$BINDING" run --help
"$PY" -I -B -S "$BINDING" prepare "$RUNTIME/no-such-provenance" "$RUNTIME/no-created-meta" --candidate 0000000000000000000000000000000000000000
test ! -e "$RUNTIME/no-created-meta"
```

The source controls passed the exact 13-name case, missing and extra names, path/version and repetition drift, physical frozen imports, and the five-case roster. Both files parsed, diagnostic CLI help exited zero, and the all-zero candidate exited one without creating metadata. This handoff covers source and in-memory controls. Actual committed-role provenance consumption and benchmark behavior still need review and root execution.

## Review round 1

The selected SDK path must equal the historical physical MacOSX27.2 SDK, and its `SDKSettings.json` pin must equal the historical path and digest. The captured R configuration must contain all six required flag names and every additional declared name. Declared values must be strings equal to the captured values. Captured empty strings are valid.

| Current source | SHA-256 |
| --- | --- |
| [Binding](namespace-bench.py) | `a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94` |
| [Source controls](namespace-bench-source-controls.py) | `eb6e7aabffeb76bf2415a2d9b8841626e5d7bd08a3a67d25e7599f588c8107d3` |
| [Prior binding snapshot](namespace-bench.review0.py) | `88b9f264c26c4e0d6c4f84a09ebf9384132c328c36358b66548bdbe6fb0cd998` |

```sh
"$PY" -I -B -S "$RUNTIME/namespace-bench-source-controls.py"
"$PY" -I -B -S -c 'import ast, sys; from pathlib import Path; [ast.parse(Path(p).read_bytes(), filename=p) for p in sys.argv[1:]]; print("PASS both Python sources parse")' "$BINDING" "$RUNTIME/namespace-bench-source-controls.py"
git apply --check "$RUNTIME/NAMESPACE-BENCH-FIX1.diff" # run from a copy of the three review-0 files
```

The source controls reproduced three invalid acceptances using the original conditions isolated from file I/O. The alternate SDK passed when both SDK pins were listed, and absent required and optional flag names passed when declared as null. After the fix, the controls and AST command both exited zero.

```text
PASS exact 13-namespace delta and rejection of missing, extra, path, version, and repeat drift
PASS authenticated physical frozen modules and preserved five-case, thirty-worker roster
PASS selected historical SDK and digest, missing/null/mismatched flags, and captured empty flag values
PASS both Python sources parse
```

The [review round 1 diff](NAMESPACE-BENCH-FIX1.diff) applies cleanly to the review-0 source set. The checks used source files and small metadata controls. Actual builds and timings remain pending in this handoff.
