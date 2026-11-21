# Prepare the first public-file diagnostic configuration

[prepare-fullcall.py](prepare-fullcall.py) is a metadata-only adapter for the frozen schema-1 runner. Its source SHA-256 is `fff62d823a35d760dfb7bed43271f9984d866bdb7726be28e38d1e9678903143` (12,669 B). It reads both real install phases, the common restoration, the observed R runtime probe, and the retained scientific smoke fixture. Static native inspection uses the Command Line Tools `otool -arch arm64 -L`. The executable portions are limited to Python metadata work and that inspection.

| Verified input | SHA-256 |
| --- | --- |
| Baseline guarded install | `f89b63cf200ebb102ddce0454bbe3b6496b2c2753ef2e90654502cb28dd69896` |
| Candidate guarded install | `26052f194c5f04712287e30484224e50f3ea77b7900c1b05082a8f1c697d253b` |
| Runtime probe | `7a8733e9b60e952379b26bd2d3a5ac235019552dccdcf56d6880bcf0cfa5054a` |
| Runtime guard receipt | `a5965333cec3d20b84fc960ee236df7f68c05abaec45bf75f45bf3d7e190e59b` |
| Smoke template | `b395193366af7a2d3c27457b22b486d05529e574b09744ed37deaebba14ff685` |

Both persisted build terminals report `completed`, raw exit 0, reaped leaders and `confirmed` closure. The outer baseline wrapper exited 1 after checking for a different success spelling. The adapter records each build's own terminal alongside its fresh archive, installed package tree and DLL hashes. Exact compiler/link commands, two normal staged load checkpoints, Apple clang 21 and SDK 27.2 markers, pinned R/Makeconf/toolchain files and the guarded R 4.5.3 runtime observation support the common provenance. The two source archive SHA-256s are `c4917d6e4276afaf8863d91754e77c00640315d61bef2528fca68bcb1b17698a` and `8b47c8164251f5380f282e92a21a1ce0686af468e83583c344bb780e95a9cd5d`.

The dependency receipt records complete per-package tree identities for the 94 common and 30 framework packages, full installed MSS role inventories, 86 arm64 native images, their direct load entries and physically readable absolute/loader-relative link file pins. The static audit records install names and system library declarations even when no ordinary on-disk file is exposed. Actual loaded namespace and DLL paths are verified in the scientific workers; retain their receipts for this step. Framework path checks use canonical paths. The reference preparation plan uses `normalizePath(.Library, mustWork = TRUE)` following the initial path assertion failure before the producer ran. [Updated protocol](../NEXT-FULLCALL.md)

Use this invocation after writing a current admission note. The private output root already holds `home`, `tmp` and the guarded runtime files. Any separately created reference directory can coexist with the new metadata. The adapter uses exclusive writes for its five files. A rerun after a partial metadata write needs a clearly accounted recovery of those specific paths.

```sh
TOP=/private/tmp/mss-fullfile-20261121
/opt/homebrew/opt/python@3.15/bin/python3.15 -I -B -S \
  "$TOP/runtime/prepare-fullcall.py" "$TOP/builds-01" "$TOP/fullcall-1m-01" \
  --admission-note "$TOP/fullcall-1m-01/admission.txt"
```

The outputs are `baseline-build.json`, `candidate-build.json`, `dependency-inventory.json`, `installations.json` and `smoke.json`. The write path checks the completed configuration using the frozen runner's schema, hash inventory, and provenance functions before publishing `smoke.json`. The subsequent scientific smoke uses `run.py run … --suite diagnostic --profile none`.

The actual inputs passed `--check` with exit 0 and the stated inventory counts. A disposable receipt with raw exit 0 and terminal status `ok` was rejected. The metadata write path and scientific smoke are reserved for the root operator.
