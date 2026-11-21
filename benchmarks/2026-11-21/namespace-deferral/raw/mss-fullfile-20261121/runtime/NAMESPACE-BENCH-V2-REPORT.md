# Namespace benchmark correction

The original diagnostic stopped after its first tiny pair. Both workers exited zero, and the outer diagnostic exited one. The retained `hard-cases-f2` TSVs contain 85 baseline namespaces and 71 candidate namespaces. Their exact baseline-only difference is the original 13 names plus `otel` 0.3.0. Both inventories agree on every shared version and all shared dependency paths. The MungeSumstats path identifies the selected installation in each role. The installed DBI load hook reached `otel` in the eager worker. The previous reference probe loaded `testthat`, whose own hook reached `otel` in both roles. Later case inventories remain unobserved. [Root-cause audit and raw evidence pins](OTEL-NAMESPACE-DELTA.md)

The [v2 entrypoint](namespace-bench-v2.py) authenticates the original binding at its physical path before executing its verified bytes. It preserves the original `__file__`, requires the exact reviewed 13-name set, and sets the explicit revised 14-name policy.

```text
AnnotationDbi DBI GenomicFeatures KEGGREST RSQLite VariantAnnotation
bit bit64 blob cachem fastmap memoise otel png
```

The existing source identity is preserved, and `namespace_policy_override` adds the actual wrapper pin, the authenticated audit pin, the prior and revised sets, the single added name and its reason. The original preparation and run identity checks seal and verify this field. The v2 entrypoint delegates both original CLI commands with the same arguments. The original parser, role/path/version checks, exact per-case repeats, first-pair stop, worker roster, timing, accounting and verification retain their implementations. [Original binding](namespace-bench.py) · [Build provenance](NAMESPACE-BUILDS-REPORT.md)

## Operator launch

Reuse the guarded normal installs of baseline `14ea7af51a19272c7792fbe582a3d42b291b4a1a` and candidate `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b`. Use actual provenance pointing to these builds and fresh benchmark admission. `META` and `RUN_PARENT` are fresh private physical roots. Launch the v2 diagnostic through the existing guard after preparation. All five cases require three new pairs and 30 fresh workers. The retained failed pair stays outside the new estimates. The full C9, 40M, warm and Windows gates remain open.

```sh
PY=/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15
RUNTIME=/private/tmp/mss-fullfile-20261121/runtime
V2="$RUNTIME/namespace-bench-v2.py"
"$PY" -I -B -S "$V2" prepare "$PROVENANCE" "$META" --candidate f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b
"$PY" -I -B -S "$V2" run "$META/diagnostic.json" "$RUN_PARENT/diagnostic" --suite diagnostic --profile none
```

## Source controls

| Artifact | SHA-256 |
| --- | --- |
| [Original binding](namespace-bench.py) | `a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94` |
| [Audit](OTEL-NAMESPACE-DELTA.md) | `1baeba297218e4d7fca22cda8706019e9093b20f9530ba9f62997e393d603e79` |
| [V2 wrapper](namespace-bench-v2.py) | `c5071a4503b583c51874127425aafb52aa2d0014f7e7fd3be33adadc2e9885b3` |
| [Source controls](namespace-bench-v2-source-controls.py) | `2b0b2c1f13ca0bb32f3ec6d265f6626f7665cfceb87a032650edb11cb4c9bdd3` |

```sh
"$PY" -I -B -S "$RUNTIME/namespace-bench-v2-source-controls.py"
"$PY" -I -B -S -c 'import ast, sys; from pathlib import Path; [ast.parse(Path(p).read_bytes(), filename=p) for p in sys.argv[1:]]; print("PASS both v2 Python sources parse")' "$V2" "$RUNTIME/namespace-bench-v2-source-controls.py"
"$PY" -I -B -S "$V2" prepare --help
"$PY" -I -B -S "$V2" run --help
```

The original policy rejected the pinned pair exactly for the extra baseline-only `otel` row. Before applying the override, the new control also failed on the baseline comparison view. All four commands above then exited zero with the revised wrapper.

```text
PASS original policy rejects actual namespace inventories (85 baseline, 71 candidate)
PASS exact 14-name policy on retained pair and rejection of omissions, extras, path/version and repeat drift
PASS original plus wrapper/audit source identity and unchanged five-case, thirty-worker roster
PASS both v2 Python sources parse
```

The controls check both raw TSV digests against the audit, use the original frozen parser, and use in-memory changes to check omissions of each expected name, extras, path/version drift and repetition drift. Both help commands printed the original argument contract. The [addition diff](NAMESPACE-BENCH-V2.diff) applied cleanly and reproduced all three added files. This is source verification. Fresh v2 preparation and measurement are pending.
