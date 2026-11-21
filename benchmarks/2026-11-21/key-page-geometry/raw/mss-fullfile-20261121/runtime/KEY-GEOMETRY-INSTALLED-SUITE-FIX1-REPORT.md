# Installed-suite post-run stop handling

The fixed operator saves the phase outcome, marks post-run verification incomplete and checks the frozen guard before reading or hashing post-run evidence. It checks for stops throughout installed-file traversal, staged-test reads, reference-link checks, runtime-control verification and receipt parsing. The existing stop-aware writer publishes `after.json`. An ordinary closed test failure can retain the verified post-run inventory while its wall budget permits. If a wall or shutdown stop interrupts verification, the funded failure result records the phase outcome and incomplete status. An inventory already published before a later stop remains available. [Source](/private/tmp/mss-fullfile-20261121/runtime/key-geometry-installed-suite.py) · [Review](/private/tmp/mss-fullfile-20261121/runtime/KEY-GEOMETRY-INSTALLED-SUITE-SOURCE-REVIEW.md)

| Artifact | SHA-256 |
| --- | --- |
| Fixed operator | `d06a1286dfcd183ea00aa1a5f2d8c73061b44bf9cbdf9a35912444aab8a75a91` |
| Original operator retained for comparison | `81dc2fb96d893137b4791521caf812dfadbe1c9ef78be6af63d0970005c3adda` |
| Review | `87ea651e681a177ee7ba21f55f3eda449c522e1c552812b8068d8e188bd1a82f` |
| Fix diff | `72e65c80e67421aff93b741d1da083fa42693c4edee89fca58856acc7aab82b5` |
| Verification controls | `00f1d63e10f5f5cf84db9c19ca59598bf9aad93b20e41b25ee9bc82cba719601` |
| Frozen guard | `58e0818aa20801c30194bdce26bd117112b4b4d5d10f74ea556360c0867b56a2` |

The controls used the production guard's actual deadline and signal check and the reused metadata writer. They substituted only the child-phase outcome and injected stops at selected evidence reads or during publication. Each control used its own temporary output directory. The actual qualification root remained absent, and syntax and single-dispatch checks passed. R qualification remains pending. [Diff](/private/tmp/mss-fullfile-20261121/runtime/KEY-GEOMETRY-INSTALLED-SUITE-FIX1.diff) · [Control script](/private/tmp/mss-fullfile-20261121/runtime/key-geometry-installed-suite-fix1-controls.py)

```sh
/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15 -I -B -S /private/tmp/mss-fullfile-20261121/runtime/key-geometry-installed-suite-fix1-controls.py
```

The command exited 0 and printed the following.

```text
Original reproduced the wall-entry defect with post-stop reads and after.json publication.
Fixed wall at entry: no post-run reads or after.json, funded failure outcome retained.
Fixed shutdown at entry: no post-run reads or after.json, funded failure outcome retained.
Shutdown during installed inventory stopped after one installed read and before after.json.
Shutdown at after.json publication was blocked by the reused stop-aware writer.
Shutdown during the remaining control checks skipped subsequent reads and marked verification incomplete.
Ordinary closed R failure retained the verified post-run inventory and a failed result.
Retained receipt census still parsed as 269 cases and 3077 passing assertions.
```
