# Key-page geometry installed-suite operator

The operator is ready for source review and one local Darwin qualification. The checks below used Python and retained evidence. The qualification command has yet to run, and its output root was absent at source verification. [Operator](/private/tmp/mss-fullfile-20261121/runtime/key-geometry-installed-suite.py)

| Input | SHA-256 |
| --- | --- |
| Suite operator | `81dc2fb96d893137b4791521caf812dfadbe1c9ef78be6af63d0970005c3adda` |
| Qualified Darwin guard | `58e0818aa20801c30194bdce26bd117112b4b4d5d10f74ea556360c0867b56a2` |
| Reused metadata writer in namespace-build producer | `615535da8374c81d3d7f23f0300c313d26d06ebea638e73239492243f851a65f` |
| Candidate build receipt | `4f4cc66733417e431281dc1892d71f5acc3fe155b2000d2404814f2c5793b74a` |
| Candidate staged install archive | `4c30721a6439428f8a7a0905a0a131304742189799c544496a0eee5e78f968ed` |
| Installed candidate package tree | `bc050421a516e28852c6e18a7fe8bacd8f7e4fbfc719327cf4e408acd20946d0` |
| Prior installed suite admission | `cbab726393fca23c7a53265bd07383d263798774848d0791b85da0034efe3ed2` |
| Prior installed suite result | `383561ec242ddbdb66d619b9ce3b826f46c783fbfb40c31de226fcf108284090` |
| Prior deferred-build guard receipt | `91c9e317497311468927b73f354b113006edb44d9a9673bb9e4dba63b8275def` |
| Prior deferred-build test-case table | `f1f878ae3c82dd9a43a538269bb76204ede85f6c70b04e76f4cdecdc8a24c333` |
| Installed-suite `control.R` | `1cc32ff00c4d8455a6d1ec7d05316f11ef0a4d6822607432effcd8c338de3c39` |
| Unfiltered installed-mode `qualify.R` | `087e55847f6a5333781f7526625ea5e59776cfd81f18f83a3c7492d55cd712b9` |
| Real union manifest | `d30550f4c01897764c858eefe03a1d3f037d3cdf0ee90be1224a419c7a2a3213` |

The candidate is the normal staged install of `8b85aa5ffc45ad7632cfc5319e85dd5eb3089be6` at `/private/tmp/mss-key-geometry-builds-20261121/candidate/library`. The build retained its stage input as an uncompressed Git archive. The operator authenticates that archive, checks its embedded commit, and records all 136 product hashes from `R/`, `src/`, DESCRIPTION and NAMESPACE. It recomputes the build's package tree and DLL hashes from the installed files. The source check measured 695 installed files and 8,220,616 bytes. [Build receipt](/private/tmp/mss-key-geometry-builds-20261121/candidate/build.json) · [Staged archive](/private/tmp/mss-key-geometry-builds-20261121/candidate/MungeSumstats.tar)

The operator copies the 92 prior staged test and support files byte for byte into the new private root, then records their hashes and whole-install hashes before and after the suite. All staged test pins match the earlier admission. Comparing the archive and staged suite found two differences, in `test-check_impute_se_beta.R` and `test-read_vcf.R`. Both change the ALS fixture path to a byte-checked copy of the installed fixture in R scratch and leave the feature calls and assertions intact. The operator checks the pinned real-reference manifest and all 102 resolved union links before and after the suite. It uses the manifest's retained package proofs and verifies the Rscript, R engine and libR pins directly. [Prior runtime review](/private/tmp/mss-fullfile-20261121/runtime/NAMESPACE-DEFERRAL-RUNTIME-REVIEW.md) · [Test admission](/private/tmp/mss-namespace-installed-suite-20261121/admission.json) · [Reference manifest](/private/tmp/mss-reference-closure-20261121/attempt-02/union-inventory.json)

The accepted guard owns the one foreground Rscript phase. The private root is `/private/tmp/mss-key-geometry-installed-suite-20261121`, with 3,600 seconds total, an 8 GiB sampled process-group RSS limit, 2 GiB of output and a continuing 186 GiB free-space floor. The operator requires floor plus cap at admission. It records the actual command, limits, root identity and minimal one-thread environment with `en_US.UTF-8`, and tracks shutdown with the accepted guard and metadata writer. The unchanged control verifies the ordered candidate, union and framework library paths. [R control](/private/tmp/mss-performance-20261119/tools/darwin-functional-stage178/control.R) · [Guard](/private/tmp/mss-performance-20261119/tools/darwin-full-stage185/guard.py)

The comparator derives 89 active files, 269 cases and 3,077 passing assertions from the prior run. Across the staged set there are 90 test scripts including one comments-only script. Python compiled the operator, checked that the source has a single guard dispatch, authenticated both reused Python controls, and recomputed the archive, install, test and reference-path checks. Receipt assertions accepted the retained successful guard record and rejected cases with nonzero R exit, unresolved closure, an unreaped leader or a missed sample. The parser recovered two existing expectation warnings and zero suite warnings from the retained run. The fresh result will record observed counts and warning totals. Completion requires a reaped R exit 0 with confirmed closure, full per-case assertion census, zero failures, skips and errors, and matching input inventories. The work here covers Darwin source readiness. Windows qualification and current R runtime results remain pending. [Prior case table](/private/tmp/mss-namespace-installed-suite-20261121/receipts/deferred/tests.tsv) · [Prior first-party qualification](/private/tmp/mss-namespace-installed-suite-20261121/receipts/deferred/qualification.tsv)

## Launch command

```sh
/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15 -I -B -S /private/tmp/mss-fullfile-20261121/runtime/key-geometry-installed-suite.py
```

Check the new root's `RESULT.json`, guard JSONL, `before.json`, `after.json` and three TSV files under `suite/` for retained run evidence.
