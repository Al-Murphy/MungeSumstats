# Retained 1M and 10M verifier handoff

The verifier is ready for source review. It checks six original editions at each scale. It reuses the original wait4 and in-worker processing timings, runs the frozen edition verifier, and starts one manifest-only R process per scale in the new evidence directory.

| Source or control | SHA256 |
| --- | --- |
| `verify-completed-scales.py` | `f6d90ac9d20af9a469f0f9098bc39ce220f8393bdb1a733cfdadd897d3edb14b` |
| `completed-scales-source.diff` | `4aa14bdab3c84432df047fd969487da3f68f9ea5d7874ef85ee23da4b1f7d2c0` |
| `completed-scales-source-controls.py` | `6e50f41fe5ebbf0397100ec1a58ada857f4bd9d7b6ee28d230a17e3677cfa73c` |
| Original cold `matrix/pins.json` | `a757eae6a2dc5f2738fb0b03f865a6bb31cb7b2a4540999bebb737e09fea6350` |
| Original failure observation | `39a96694ac60f09bfad5e662c24ec2759fc88823149b2e94d2e3bd1628ecf995` |
| Original configuration | `156143e169744ead89eb7a01721ac433da0d0983be68924b417bcc72f00426c7` |
| Original SNP profile | `b387920dbab4d53af316a8b2b320d577be463abdad7ce80dbf8b9de411e303d0` |
| Authenticated `warm-current.py` | `c58b14cda2ad650cefa078e4a73cd2a2a4d5ab8ffc4f7477e5375e0aaa794220` |

The full addition diff is `/private/tmp/mss-fullfile-20261121/runtime/completed-scales-source.diff`.

## Source-only results

```sh
py=/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15
runtime=/private/tmp/mss-fullfile-20261121/runtime
"$py" -I -B -S "$runtime/completed-scales-source-controls.py"
"$py" -I -B -S "$runtime/verify-completed-scales.py" --help
/usr/bin/shasum -a 256 "$runtime/verify-completed-scales.py" "$runtime/completed-scales-source-controls.py" "$runtime/completed-scales-source.diff"
git diff --no-index -- /dev/null "$runtime/verify-completed-scales.py" > "$runtime/completed-scales-source.diff"
git diff --no-index --check -- /dev/null "$runtime/verify-completed-scales.py"
```

The controls returned 0. They compiled the verifier, authenticated the helper imports, selected the two rows and pair counts, rejected a wrong pair count and a false original terminal record, and rejected a missing worker directory and a wrong configuration hash. `--help` returned 0. The diff command returned the expected change status 1 and the whitespace check emitted no diagnostics. Source and small control reads were the only checks. Scientific execution, R, payload hashing, and installed-tree traversal remain for the guarded run.

## Guarded invocation for the operator

Create `/private/tmp/mss-completed-scales-20261121/verify-01` as a physical, operator-owned directory with mode 0700. Use the reviewed source pin, existing sanitized one-thread environment and diagnostic guard, and dispatch this payload inside that guard with allowed exit 0.

```sh
/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15 \
  -I -B -S /private/tmp/mss-fullfile-20261121/runtime/verify-completed-scales.py \
  --output /private/tmp/mss-completed-scales-20261121/verify-01
```

The scope uses the proposed 4 GiB RSS, 16 MiB own-root metadata, 186 GiB filesystem floor and 600 s wall envelope. The original opening records about 2.66 GB of scientific assets per inventory pass. The two inventory passes, installed tree checks, and the retained edition scans are untimed verifier work. Their total duration has yet to be measured.

## Expected evidence

- `current-verification-pins.json` records fresh UTC, the original opening and the live asset, toolchain, role and dependency inventories.
- `matched-1M/` and `matched-10M/` each contain digest records for every artifact in six retained editions and the stable RDS manifest result and process controls.
- `scoped-summary.json` retains per-pair ratios, processing medians, whole-worker wait4 wall and CPU medians, RSS median and maximum, relative reductions, and the separate 8 GiB candidate evaluation. Exit 0 indicates successful retained-edition verification with the performance status recorded separately.
- `oracle-window-plan.json` keys the windows by `1M` and `10M`. Each window binds the original pair-01 candidate edition by physical path and device/inode, the `mss7=` selector, the artifact digest control, and full SHA256/MD5/byte digests for the input and three edition files.
- `RESULT.json` pins the verifier source, helper sources, original small controls, selected verification evidence, and the original disk-floor failure with raw exit -15 and supervisor exit 1. It records the final equality of the current inventories and flags the two independent oracles and root observer as pending.

Each independent 19-field oracle and its before/after digest window require separate guarded execution and observer closure. The original full matrix remains incomplete and warm ineligible. The retained 40M baseline, censored candidate, supplemental cases and C9 are outside this result.

## Review round 1 fix

The worker collector now uses authenticated `warm.ROLE_LOGS` to pin stdout, stderr and `r-heap.rds` for each of the twelve selected workers. The original review source is retained as `verify-completed-scales.review0.py`. The top source table records that initial review state. Use the current source pin below for dispatch.

| Source or control | SHA256 |
| --- | --- |
| Current `verify-completed-scales.py` | `fb61736c23dcde4bfe16fedfa450c72a67d8803273536022193d4d39004e2848` |
| `verify-completed-scales.review0.py` | `f6d90ac9d20af9a469f0f9098bc39ce220f8393bdb1a733cfdadd897d3edb14b` |
| `COMPLETED-SCALES-FIX1.diff` | `12335cb92cc6640992b436e9c9283ba315e9461ce488509f1963137a1631dfcf` |
| Current `completed-scales-source-controls.py` | `0b35f90be7285ed916057193d2d7d09f2297e778c7b1fb0793af8bc3c64ad928` |

```sh
py=/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15
runtime=/private/tmp/mss-fullfile-20261121/runtime
"$py" -I -B -S "$runtime/completed-scales-source-controls.py"
"$py" -I -B -S "$runtime/verify-completed-scales.py" --help
git diff --no-index -- "$runtime/verify-completed-scales.review0.py" "$runtime/verify-completed-scales.py" > "$runtime/COMPLETED-SCALES-FIX1.diff"
git diff --no-index --check -- "$runtime/verify-completed-scales.review0.py" "$runtime/verify-completed-scales.py"
/usr/bin/shasum -a 256 "$runtime/verify-completed-scales.review0.py" "$runtime/verify-completed-scales.py" "$runtime/completed-scales-source-controls.py" "$runtime/COMPLETED-SCALES-FIX1.diff"
```

The new in-memory collector control exited 1 before the one-line fix and identified the missing heap diagnostic.

```text
syntax: pass
authenticated source imports: pass
in-memory roster selection and wrong pair count rejection: pass
failure binding and false terminal record rejection: pass
missing worker rejection: pass
wrong configuration digest rejection: pass
Traceback (most recent call last):
  File "/private/tmp/mss-fullfile-20261121/runtime/completed-scales-source-controls.py", line 64, in <module>
    assert not missing, f'Missing worker log pins: {sorted(map(str, missing))}'
           ^^^^^^^^^^^
AssertionError: Missing worker log pins: ['/synthetic/pair-01-baseline/r-heap.rds']
```

The same source-only command exited 0 after the fix.

```text
syntax: pass
authenticated source imports: pass
in-memory roster selection and wrong pair count rejection: pass
failure binding and false terminal record rejection: pass
missing worker rejection: pass
wrong configuration digest rejection: pass
worker log pins include the heap diagnostic: pass
```

`--help` exited 0. The review diff contains one changed line, and the whitespace check emitted no diagnostics. This round used in-memory worker timing and log-pin records. Guarded scientific work remains pending.
