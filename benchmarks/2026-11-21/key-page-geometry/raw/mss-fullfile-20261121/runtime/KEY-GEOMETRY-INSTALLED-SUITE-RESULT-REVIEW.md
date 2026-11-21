# Darwin installed-suite result review

**Specification compliance — pass. Evidence quality — pass.** The normal staged installation of candidate `8b85aa5ffc45ad7632cfc5319e85dd5eb3089be6` passed the unfiltered installed testthat suite on Darwin arm64. Open findings are 0 Critical, 0 Important and 0 Minor. Native Windows qualification and performance measurements are outside this result.

| Observed from the test receipts | Result |
| --- | ---: |
| Active test scripts | 89 of 90 |
| Test cases | 269 |
| Passing assertions | 3,077 |
| Skips, failures and errors | 0 each |
| Test warnings | 2 |
| Suite warnings | 0 |

I recounted the case table and byte-compared all three TSV receipts with those from the prior deferred-build suite. The case table, conditions and qualification summary matched exactly. The remaining test script contains only comments. The two warnings occurred in the same cases as before, for a replaced import in the a0/a1 test and numeric coercion in the rs-ID test. [Cases](/private/tmp/mss-key-geometry-installed-suite-20261121/suite/tests.tsv) · [Conditions](/private/tmp/mss-key-geometry-installed-suite-20261121/suite/conditions.tsv) · [Qualification](/private/tmp/mss-key-geometry-installed-suite-20261121/suite/qualification.tsv)

The guard's planned command matches the admission and invokes the frozen suite control on the candidate library with the 102-package reference union. Admission identifies the reviewed operator, candidate commit, archived source tree and 136 product paths. The before and after inventory records are byte-identical, cover 695 installed files and 92 staged sources, and bind the installed package tree to the build receipt. The staged source digests also match the prior admission. The operator recorded complete post-run verification. [Admission](/private/tmp/mss-key-geometry-installed-suite-20261121/admission.json) · [Before](/private/tmp/mss-key-geometry-installed-suite-20261121/before.json) · [After](/private/tmp/mss-key-geometry-installed-suite-20261121/after.json)

The guard recorded a reaped raw exit of 0 with confirmed closure, zero missed samples, no signals and no cleanup errors. The phase lasted 170.208 seconds against the 3,600-second window. Peak sampled group RSS was 6,612,402,176 bytes against 8 GiB. Peak allocated output was 19,185,664 bytes against 2 GiB, and the minimum observed free space was 230,245,568,512 bytes against the 186 GiB floor. The independent observation records outer exit 0 and no surviving member of group or session 78054, with zero inspection errors. The guard terminal, operator result and observation agree. [Guard transcript](/private/tmp/mss-key-geometry-installed-suite-20261121/candidate-unfiltered.jsonl) · [Operator result](/private/tmp/mss-key-geometry-installed-suite-20261121/RESULT.json) · [Independent observation](/private/tmp/mss-key-geometry-installed-suite-observation-20261121.json)

| Evidence | Verified SHA-256 |
| --- | --- |
| Operator | `d06a1286dfcd183ea00aa1a5f2d8c73061b44bf9cbdf9a35912444aab8a75a91` |
| Operator result | `c1d6da57e60e552b3b17c6f07e88363989962ce8b5c074450f77498b58c7770c` |
| Before and after inventory records | `33304938699806fc08c223599b62802f2a332c2cd6fca016024302dc9e770505` |
| Independent observation | `4b099889a6a39c13d0b2950bcbbc14d0d0f868113dc789cfc3c50f508bde7017` |

This review used the retained records from the actual run. Installed-file and archive measurements come from the reviewed operator. I verified their recorded identities and the small evidence files without rerunning R or rehashing package payloads.
