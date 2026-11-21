# Key-page geometry public-API benchmark binding

The 16 KiB key-page trial merits a matched public-API comparison after preserving native results and reducing nominal requested key bytes at the same 8 MiB key payload budget. The trial increased application seek/read calls. Paired public-API elapsed time and peak RSS remain unmeasured. [Independent native counter review](KEY-GEOMETRY-RESULT-REVIEW.md)

The [wrapper](key-geometry-bench.py) authenticates the original namespace binding at its physical path before executing those verified bytes. It checks the original `14ea7af51a19272c7792fbe582a3d42b291b4a1a` baseline and exact prior 13-name namespace policy, then sets the deferred baseline to `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b` and the expected baseline-only set to empty. Both roles defer `VariantAnnotation`. The original callback therefore requires both roles to load identical namespaces within each case, with matching versions and physical dependency paths and exact repetition within each role. The original parser still verifies each role's selected MungeSumstats installation before normalizing its role-specific path. A first-pair mismatch stops the run before pair two. The retained 71-name inventory is the source control. [Original binding](namespace-bench.py) · [Retained deferred inventory](/private/tmp/mss-namespace-benchmark-v2-20261121/diagnostic/hard-cases-f2/pair-01-candidate/namespaces.tsv)

The wrapper preserves the original binding's physical `__file__` and source identity, with the new baseline reflected by the original `baseline_override` field. It adds `key_geometry_override` with the actual wrapper pin, original and actual baseline IDs, prior and revised namespace policies, the reason for the policy, and the checked independent counter-review pin. The original preparation and run identity checks seal and verify this field. The original and v2 namespace sources retain their existing hashes.

Low-level build and runtime provenance, pinned dependencies and assets, the scientific case selection, worker timing and `wait4` accounting, tiny literal golds, stable edition bytes and RDS parity retain their original implementations. All five expanded cases use three pairs, giving 30 fresh workers in unprofiled diagnostic mode. The cases are `hard-cases-f2`, `hard-cases-f3`, `declared-f2`, `declared-f3`, and `matched-1M`.

## Operator inputs and launch

Supply fresh matched normal builds for the deferred baseline and the actual committed trial, with complete original provenance. Preparation requires the distinct actual 40-hex candidate commit. Use fresh private roots for preparation and measurement. Admit the planned 3 GiB retained output and 186 GiB remaining free floor using fresh host observations. The geometry trial uses its own cohort. Run the independent 1M oracle over all 19 fields after timing. Full C9, full matrix, warm, native Windows and publication remain open. [Original receipt contract](NAMESPACE-BENCH-REPORT.md)

```sh
PY=/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15
RUNTIME=/private/tmp/mss-fullfile-20261121/runtime
BINDING="$RUNTIME/key-geometry-bench.py"
"$PY" -I -B -S "$BINDING" prepare "$PROVENANCE" "$META" --candidate "$ACTUAL_CANDIDATE_COMMIT"
"$PY" -I -B -S "$BINDING" run "$META/diagnostic.json" "$RUN_PARENT/diagnostic" --suite diagnostic --profile none
```

## Source verification

| Artifact | SHA-256 |
| --- | --- |
| [Wrapper](key-geometry-bench.py) | `c2eb2af6b254a981272740eb05b6632b4dc0b55603f62c6cb9a1149a62dee976` |
| [Source controls](key-geometry-bench-source-controls.py) | `4d9913d3424514becea1a066187db75ce59e385195a6f9fbb59364435cde23d0` |
| [Original binding](namespace-bench.py) | `a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94` |
| [Namespace v2 wrapper](namespace-bench-v2.py) | `c5071a4503b583c51874127425aafb52aa2d0014f7e7fd3be33adadc2e9885b3` |
| [Native counter review](KEY-GEOMETRY-RESULT-REVIEW.md) | `d478ab3a248afb1f94cca78c913037735aadea29d8e6bd9c4fbbd539c4f425d0` |
| [Retained f39 worker TSV](/private/tmp/mss-namespace-benchmark-v2-20261121/diagnostic/hard-cases-f2/pair-01-candidate/namespaces.tsv) | `e6c40018fe19afb387fdb4fe2c0bc2a4f30302208fca91eb56364ac943980c32` |

```sh
"$PY" -I -B -S "$RUNTIME/key-geometry-bench-source-controls.py"
"$PY" -I -B -S -c 'import ast, sys; from pathlib import Path; [ast.parse(Path(p).read_bytes(), filename=p) for p in sys.argv[1:]]; print("PASS both key-geometry Python sources parse")' "$BINDING" "$RUNTIME/key-geometry-bench-source-controls.py"
"$PY" -I -B -S "$BINDING" prepare --help
"$PY" -I -B -S "$BINDING" run --help
```

The control read the pinned TSV through the original frozen parser. The original eager gate rejected it. The deferred policy accepted two equal copies of that normalized state and rejected the in-memory perturbations. All four commands exited zero with the completed wrapper.

```text
PASS the retained 71-name inventory matches itself and missing, extra, version, path and repeat drift fail
PASS original, wrapper and counter review pins, deferred baseline, and unchanged five-case, thirty-worker roster
PASS both key-geometry Python sources parse
```

Both help commands printed the original argument contract. The [addition diff](KEY-GEOMETRY-BENCH.diff) applied cleanly and reproduced all three added files. Fresh geometry builds, preparation, measurement and the independent oracle are pending.
