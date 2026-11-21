# Matched committed namespace builds

The source is ready for independent review and execution. Both role installations and all benchmark observations are still pending.

| Role | Commit |
| --- | --- |
| Baseline | `14ea7af51a19272c7792fbe582a3d42b291b4a1a` |
| Candidate | `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b` |

The producer takes a single absent absolute output path under a physical parent. It creates its own private root and records fresh filesystem availability in `admission.json`. It admits the independent build when the initial free space funds the 186 GiB continuing floor and the 512 MiB owned output budget. One authenticated native guard shares the 900-second operator window across all phases, samples group RSS against 4 GiB, and checks logical and allocated output against the budget. The operator uses the native guard's exclusive files and records its actual returned outcome alongside each retained terminal. [Producer](namespace-builds.py) · [Guard](/private/tmp/mss-performance-20261119/tools/darwin-full-stage185/guard.py)

The producer authenticates the reviewed binding before executing its bytes and obtains its frozen helpers from `load_frozen()`. It checks the original matrix and dependency inventory pins, freshly verifies their compiler, SDK, R Makeconf and runtime file pins, and captures the pinned clang's `--version` plus effective `R CMD config --all` under the guard. The exact frozen worker captures `runtime.tsv` with no selected package library. The actual runtime must match the original matrix, including R, BLAS, LAPACK and locale. The toolchain object is shared by both receipts and retains the captured compiler lines and required flags, including any explicit empty values. [Binding](namespace-bench.py) · [Provenance consumer contract](NAMESPACE-BENCH-REPORT.md)

Each role runs a guarded local Git tree lookup and archive from `/private/tmp/mss-lazy-variantannotation-20261121`. The guard owns the archive's stdout file `ROLE/MungeSumstats.tar`, directly receiving the bytes from the requested `git archive --format=tar --prefix=MungeSumstats/ COMMIT`. The `/usr/bin/env` wrapper applies the explicit Git safety settings without modifying the authenticated guard. The operator result retains the wrapper dispatch, direct command, effective environment, raw outcome, transcript pins and actual terminal. The benchmark consumer independently checks the embedded archive commit.

The two installs use the same clean build environment with private HOME and temporary paths, the physical 96-package common library, disabled user startup and Makevars, one make job and the guard's supported thread controls. The approved `en_US.UTF-8` locale applies to both installs and to the runtime capture. It matches the frozen worker's expected runtime TSV and the binding's exact SDK log marker. Both installs run the same SDK-prefixed `R CMD INSTALL --use-vanilla` invocation with separate libraries and archives, without an outer timing wrapper. The producer checks staged load markers and three parsed native commands, then records DLL and frozen package-tree digests. It requires distinct source trees and identical native commands before publishing `PROVENANCE.json`. [Build requirements](NAMESPACE-BUILDS-BRIEF.md)

The provenance pins each role's `build.json` and the fresh runtime TSV, with the runtime's guard JSONL, stdout and stderr in that order. Both role receipts use the same transcript order. Provenance declares the 186 GiB continuing floor and the 3 GiB estimate requested by the consumer. Benchmark admission, host contention and cache observations, and timings require their separate measurement phase. The private `operator-result.json` records only attempted phases. A failure terminates the sequence and preserves native guard evidence. The producer uses the same guard funding and root checks for its own JSON controls. It reports a write refusal to stderr if the guard cannot fund the failure result.

## Operator command

```sh
PY=/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15
RUNTIME=/private/tmp/mss-fullfile-20261121/runtime
"$PY" -I -B -S "$RUNTIME/namespace-builds.py" "$OUTPUT"
```

After a zero exit, a completed operator result and independent observation, `$OUTPUT/PROVENANCE.json` is the input to the separately reviewed binding's preparation command. [Preparation and measurement commands](NAMESPACE-BENCH-REPORT.md)

## Initial source checks

| Source | SHA-256 |
| --- | --- |
| [Producer](namespace-builds.py) | `41da5092b5aafc65c7d95cbd684e53f4606e131b5f0a818293548080d502f4df` |
| [Metadata controls](namespace-builds-source-controls.py) | `56dbbb35d66cb4481f583a17120e59207c60426833ac8f15b64c95c1167b12c0` |

The physical Python parsed both Python files and printed CLI help with exit zero. In-memory controls passed for explicit empty compiler flags, missing required flags, malformed and duplicate runtime keys, runtime locale mismatches and guard environment compatibility. Passing the existing runtime directory as output exited one at exclusive root creation. Source verification used metadata only.

## Review round 1

Every operator JSON write now applies the frozen guard's current floor, owned budget, evidence reservation and directory checks. A rejected admission with available space for its record retains `admission.json` with `admitted=false` and exits one. A rejection at the continuing floor reports the funding refusal to stderr. After a failed phase, a funded failure result retains the attempted phase records. On refusal, stderr reports the reason and prior evidence remains in place. [Review](NAMESPACE-BUILDS-REVIEW.md)

Provenance and final completion remain inside the guard's shared deadline and signal scope. The writer rechecks the guard after the directory walk and after closing the new control. A stop during the write removes the control it just created after confirming its inode and exits nonzero through the failure path. A stop after a valid provenance receipt can leave that receipt alongside a failed or absent operator result. The zero exit and completed operator result remain the handoff criteria.

| Current source | SHA-256 |
| --- | --- |
| [Producer](namespace-builds.py) | `615535da8374c81d3d7f23f0300c313d26d06ebea638e73239492243f851a65f` |
| [Metadata controls](namespace-builds-source-controls.py) | `b757ed4a868217c3476a4b8c824671cb2e35116d7942e63f9d2193ef6e606b99` |
| [Review snapshot](namespace-builds.review0.py) | `41da5092b5aafc65c7d95cbd684e53f4606e131b5f0a818293548080d502f4df` |

```sh
"$PY" -I -B -S "$RUNTIME/namespace-builds-source-controls.py"
"$PY" -I -B -S -c 'import ast, sys; from pathlib import Path; [ast.parse(Path(p).read_bytes(), filename=p) for p in sys.argv[1:]]; print("PASS both producer Python sources parse")' "$RUNTIME/namespace-builds.py" "$RUNTIME/namespace-builds-source-controls.py"
"$PY" -I -B -S "$RUNTIME/namespace-builds.py" --help
"$PY" -I -B -S "$RUNTIME/namespace-builds.py" "$RUNTIME"
git -C "$RUNTIME" apply --reverse --check NAMESPACE-BUILDS-FIX1.diff
```

```text
PASS required and empty compiler flags, strict runtime TSV, and native guard environment compatibility
PASS guarded failure and rejection funding, deadline after walk, and stop during publication
PASS both producer Python sources parse
[Errno 17] File exists: '/private/tmp/mss-fullfile-20261121/runtime'
```

The existing-root probe exited one. The metadata, syntax, help and reverse diff checks exited zero. The funding controls used the native guard's funding arithmetic and small temporary controls in the runtime directory. They exercised the owned cap, the free floor, funded and unfunded initial rejection records, deadline expiry during a walk and shutdown during publication. No guard child phase or native build ran.
