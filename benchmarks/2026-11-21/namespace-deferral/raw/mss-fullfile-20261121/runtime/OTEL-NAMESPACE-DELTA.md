# OpenTelemetry namespace delta

Both first-pair `hard-cases-f2` workers exited 0. The diagnostic supervisor stopped when the strict namespace policy found `otel` 0.3.0 as a fourteenth baseline-only namespace. The 13-name expectation came from a reference probe that loaded `testthat` before collecting its namespace inventory. The installed DBI and `testthat` load hooks independently look up `otel`. In the benchmark worker, only the eager baseline reaches DBI. In the reference probe, `testthat` supplies the second path in both roles. [Diagnostic error](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/failure.json) · [Observed worker outcomes](/private/tmp/mss-namespace-benchmark-failure-observation-20261121.json) · [13-name policy](/private/tmp/mss-fullfile-20261121/runtime/namespace-bench.py:30)

## Installed dependency path

- The installed baseline `MungeSumstats` imports eleven symbols from `VariantAnnotation`. All eleven imports are absent from the candidate `NAMESPACE`. Both `DESCRIPTION` files retain the `VariantAnnotation` declaration in `Imports`. [Baseline NAMESPACE](/private/tmp/mss-namespace-builds-20261121/baseline/library/MungeSumstats/NAMESPACE:49) · [Candidate NAMESPACE](/private/tmp/mss-namespace-builds-20261121/candidate/library/MungeSumstats/NAMESPACE:49) · [Candidate DESCRIPTION](/private/tmp/mss-namespace-builds-20261121/candidate/library/MungeSumstats/DESCRIPTION:34)
- `VariantAnnotation` imports directly from `DBI`, `AnnotationDbi`, and `GenomicFeatures`. `GenomicFeatures` also imports `DBI` and `AnnotationDbi`. `AnnotationDbi` imports `DBI`, `RSQLite`, and `KEGGREST`. `RSQLite` reaches `bit64`, `blob`, `DBI`, and `memoise`. `bit64` imports `bit`, `memoise` imports `cachem`, `cachem` imports `fastmap`, and `KEGGREST` imports `png`. These edges account for the 13 statically connected baseline-only packages in the original policy. [VariantAnnotation NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/VariantAnnotation/NAMESPACE:23) · [GenomicFeatures NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/GenomicFeatures/NAMESPACE:5) · [AnnotationDbi NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/AnnotationDbi/NAMESPACE:4) · [RSQLite DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/RSQLite/DESCRIPTION:37) · [bit64 NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/bit64/NAMESPACE:299) · [memoise NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/memoise/NAMESPACE:16) · [cachem NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/cachem/NAMESPACE:10) · [KEGGREST NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/KEGGREST/NAMESPACE:4)
- Across the installed common96, `DESCRIPTION` lists `otel` in `Suggests` for `DBI` 1.3.0 and `testthat` 3.3.2, and in `Imports` for `callr` 3.8.0. The common96 `NAMESPACE` files have zero explicit `otel` imports. The installed DBI and `testthat` `.onLoad` bodies each call their own `otel_cache_tracer()`, whose first guard includes `requireNamespace("otel", quietly = TRUE)`. That DBI hook explains the worker delta. The `testthat` hook accounts for `otel` in both probe roles. `testthat` declares `callr` in `Imports`, while all six benchmark and probe snapshots have zero `callr` and zero `otelsdk` entries. [DBI DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/DBI/DESCRIPTION:19) · [DBI release notes](/private/tmp/mss-parallel-support-20261121/common/DBI/NEWS.md:7) · [testthat DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/testthat/DESCRIPTION:18) · [testthat release notes](/private/tmp/mss-parallel-support-20261121/common/testthat/NEWS.md:3) · [callr DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/callr/DESCRIPTION:19)

### Installed hook evidence

The installed compact lazy-load index `common/<package>/R/<package>.rdx` maps each binding to a byte range in `common/<package>/R/<package>.rdb`. Read the gzip index, with a 32 KiB output cap, and decode the big-endian XDR vector of `(offset, length)` pairs and the corresponding names. For each row below, seek to the listed offset, read the listed number of bytes, and inflate the zlib payload after the four-byte uncompressed-length prefix. Each cap applies to the selected block. This procedure inspects the serialized source expressions embedded in the installed bytecode using Python's standard library. It requires no R execution.

| Package | Binding | Byte offset | Bytes read | Bytes inflated | SHA-256 of bytes read |
|---|---|---:|---:|---:|---|
| DBI | `.onLoad` | 300957 | 207 | 372 | `0f9ed4a1cf4eb637d98ec790a27af2df4d4795c9d20f696cecc87ab7197af84e` |
| DBI | `otel_cache_tracer` | 384113 | 474 | 1605 | `c9a33c5442346537a8a9c1e23571bb2cb82bd6804490b3c0933fa00af9fc647b` |
| testthat | `.onLoad` | 8326 | 209 | 377 | `da3442f140d6152537659b551285d20fc83f821d32d70e296ad3946dbf7cbfe1` |
| testthat | `otel_cache_tracer` | 1109922 | 473 | 1604 | `651801d227917809de11cbf15480aeb9e720afaa20875710f77e85016ffdc63e` |

In the inflated `.onLoad` blocks, the call head `otel_cache_tracer` starts at byte 230 in DBI and 235 in `testthat`. It is followed by the nil argument list in each. In the inflated tracer blocks, the call head `requireNamespace` starts at byte 362 in DBI and 361 in `testthat`. Each call has a length-one character vector containing `otel` and a second argument tagged `quietly` followed by an XDR logical vector of length one and value one (`0000000a 00000001 00000001`). Both blocks then contain a `return` expression and the tracer calls `otel::get_tracer` and `otel::tracer_enabled`. The independently pinned installed hook expressions and loaded-namespace observations support the two load paths above.

## Inventory comparison

The historical probe executes `loadNamespace("MungeSumstats")`, checks whether `VariantAnnotation` loaded, executes `loadNamespace("testthat")`, then records its `before` namespace inventory. `after` adds exactly six reference namespaces in each role, with no removals or changes to existing records. The benchmark worker attaches `Rcpp` before `MungeSumstats` and writes the namespace inventory after processing. Both workflows use the installed common96 and R 4.5.3. [Reference probe load sequence](/private/tmp/mss-reference-closure-20261121/attempt-03/namespaces.R:57) · [Benchmark worker load sequence](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R:100) · [Benchmark snapshot](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R:222)

| Snapshot | Eager or baseline | Deferred or candidate | Shared names | Eager or baseline only | Deferred or candidate only |
|---|---:|---:|---:|---:|---:|
| Diagnostic worker, after processing | 85 | 71 | 71 | 14 | 0 |
| Reference probe, before reference loads | 87 | 74 | 74 | 13 | 0 |
| Reference probe, after reference loads | 93 | 80 | 80 | 13 | 0 |

The raw worker rows agree exactly on versions for all 71 shared names and on physical paths for 70. The sole path difference is the selected `MungeSumstats` 1.21.4 installation under `namespace-builds-20261121/{baseline,candidate}/library/MungeSumstats`. The 85 baseline rows comprise 69 common packages, 15 framework packages, and the role installation. The 71 candidate rows comprise 55 common packages, the same 15 framework packages, and the role installation. All 14 baseline-only packages resolve physically under `/private/tmp/mss-parallel-support-20261121/common/<package>`.

| Worker baseline only | Raw version |
|---|---|
| AnnotationDbi | 1.72.0 |
| DBI | 1.3.0 |
| GenomicFeatures | 1.62.0 |
| KEGGREST | 1.50.0 |
| RSQLite | 3.53.3 |
| VariantAnnotation | 1.56.0 |
| bit | 4.6.0 |
| bit64 | 4.8.6 |
| blob | 1.3.0 |
| cachem | 1.1.0 |
| fastmap | 1.2.0 |
| memoise | 2.0.1 |
| otel | 0.3.0 |
| png | 0.1.9 |

The 13-name reference-probe difference consists of those rows except `otel`. Both probe stages record `otel` 0.3.0 at `/private/tmp/mss-parallel-support-20261121/common/otel` in both roles. The reference-probe eager `before` adds just `brio` 1.1.5 and `testthat` 3.3.2 relative to the benchmark baseline. The reference-probe deferred `before` adds just `brio`, `testthat`, and `otel` relative to the benchmark candidate. Both sets use the same common physical paths. The six `after` additions are `BSgenome.Hsapiens.1000genomes.hs37d5` 0.99.1, `BSgenome.Hsapiens.NCBI.GRCh38` 1.3.1000, both `SNPlocs.Hsapiens.dbSNP144` builds at 0.99.20, and both `SNPlocs.Hsapiens.dbSNP155` builds at 0.99.24, each physically identical between the probe roles.

For cross-workflow comparisons of matching roles, the only differences in the shared rows are the selected `MungeSumstats` paths and version punctuation. Probe eager uses `/private/tmp/mss-namespace-annfix-build-20261121/library/MungeSumstats`, and probe deferred uses `/private/tmp/mss-namespace-deferred-build-20261121/library/MungeSumstats`. `getNamespaceVersion()` in the probe preserves hyphens that `packageVersion()` in the benchmark renders as dots. The raw spellings differ for `Matrix` (`1.7-4`), `RCurl` (`1.98-1.20`), `XML` (`3.99-0.25`), `abind` (`1.4-8`), `bitops` (`1.1-0`), `codetools` (`0.2-20`), `lattice` (`0.22-9`), and eager-only `png` (`0.1-9`). After normalizing `-` to `.`, all shared versions agree across workflows. Within the reference probe, every shared version and all shared physical dependency paths agree across its roles at each stage.

The exact 14-name difference is established for this completed pair. The first-pair stop leaves later case inventories unobserved in this diagnostic.

## Raw evidence pins

SHA-256 digests cover the full small files in this table. The `.rdb` hashes above cover only the specified byte ranges.

| Evidence file | Bytes | SHA-256 |
|---|---:|---|
| [Diagnostic error](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/failure.json) | 122 | `1c57528fb03fdcfeb9caede718598fbcc28cb1165004e19924d2da3b30512737` |
| [Diagnostic input pins](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/pins.json) | 402806 | `bca85f72b4763fa93a462b3109402179b1f64d388c93d1e85cee0851b754ac72` |
| [Worker baseline namespaces](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/hard-cases-f2/pair-01-baseline/namespaces.tsv) | 7111 | `040814dc3654a620b7096ab118e8e72d2d75a9204a77cc849f96b6b501fff255` |
| [Worker baseline process](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/hard-cases-f2/pair-01-baseline/process.json) | 713 | `5b1f60c5ba67b7194fa0f124fc31c56e62064d9b5dcbe1e5200c2862abfef45e` |
| [Worker candidate namespaces](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/hard-cases-f2/pair-01-candidate/namespaces.tsv) | 6007 | `e6c40018fe19afb387fdb4fe2c0bc2a4f30302208fca91eb56364ac943980c32` |
| [Worker candidate process](/private/tmp/mss-namespace-benchmark-20261121/diagnostic/hard-cases-f2/pair-01-candidate/process.json) | 715 | `b51ca1e4fa95001776d7c68cfafdc3f7e98e348fdf55d11bb604f2363a1ea4d0` |
| [Reference eager namespaces](/private/tmp/mss-reference-closure-20261121/attempt-03/eager-namespaces.json) | 52522 | `f366b1a7083252c2593baa0da5b525fcc157171e3f21225c4e8e221c8ce80174` |
| [Reference deferred namespaces](/private/tmp/mss-reference-closure-20261121/attempt-03/deferred-namespaces.json) | 45632 | `e452937dc8885bf9175db4ab32af4feb68d074624accdbdcfc3b1785bf8c4fd6` |
| [Reference probe](/private/tmp/mss-reference-closure-20261121/attempt-03/namespaces.R) | 4282 | `279ab8869dc2d35c68a4d9d79fcc8dc7b81ffe087136f5989d2b02b708e89e14` |
| [Benchmark worker](/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade/worker.R) | 12721 | `f8477b5b0047b0361d636e4af63b5ed9bf2be87229d7cf9c682331a51147e7df` |
| [Benchmark namespace policy](/private/tmp/mss-fullfile-20261121/runtime/namespace-bench.py) | 24142 | `a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94` |
| [DBI DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/DBI/DESCRIPTION) | 1947 | `0b1bde3978d5641d47d6ca86870e9feb6a8d27d29aff82f100c88c6167c4fc5e` |
| [DBI NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/DBI/NAMESPACE) | 3280 | `7be0aab013580a78b030d2dc51759a3dabfaaa7ffa1f100bec3cfd853e8828b0` |
| [DBI lazy-load index](/private/tmp/mss-parallel-support-20261121/common/DBI/R/DBI.rdx) | 6142 | `768383cbf4e167266f086fbe840ff2ac9fe7fa6d6ff76f3be0fdd0c89800d6f7` |
| [testthat DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/testthat/DESCRIPTION) | 1877 | `67c7d804cc21f8b29691d55e9a8d7120691bb07f348ed9d575e9fa5d7891d5d3` |
| [testthat NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/testthat/NAMESPACE) | 5997 | `88f3b377dc878144b6278db57d76071ee8a8089e37cf5b9540af253a11c81de5` |
| [testthat lazy-load index](/private/tmp/mss-parallel-support-20261121/common/testthat/R/testthat.rdx) | 8493 | `ae21972442842b3734900e1aeecf2a4f62af04be2a3bd9aafd27b7e08eef3d9d` |
| [callr DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/callr/DESCRIPTION) | 1657 | `c02a5a8def43dc911a8c428c2527f0dac9316f4ffa0f2b60f88f8c6d60365dbe` |
| [callr NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/callr/NAMESPACE) | 604 | `88d441671c169492928077bfbe6a5bfd97f0f2df558ca2ba1a3d3b57c978f43f` |
| [otel DESCRIPTION](/private/tmp/mss-parallel-support-20261121/common/otel/DESCRIPTION) | 1357 | `c8eee53dbdca80909aee62fff04309049cf825e06af18d88d7866f37ec1e091b` |
| [otel NAMESPACE](/private/tmp/mss-parallel-support-20261121/common/otel/NAMESPACE) | 1758 | `aaba98588b9f20339b1914e48b1f1f30ef5a6e940940de78f9d1ad861764341a` |
