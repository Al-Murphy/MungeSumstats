# Complete installed-package CI

The **Installed package qualification** workflow builds the commit that triggered it, installs the generated source package with ordinary R, and tests that installation. The existing fast standalone workflow remains available. The two workflows test different scopes.

The full producer uses standard hosted Ubuntu22.04 and Windows2022 x64 runners, Python3.12.10, R4.5.3, Bioconductor3.22 and Rtools45 on Windows. Every action is pinned to an official repository commit. Repository/run/attempt/job identities come from GitHub, rather than a maintainer account or private carrier branch.

## What runs

1. Linux normally builds and installs the triggering checkout and generates genuine outputs for allele flipping, indel retention, multiple rsIDs, ordered/reversed tabix queries and plain/gzip/BGZF TSV/VCF headers. It records the actual installed source/DLL, complete outputs and R column types. It runs the existing Linux pidfd/subreaper ownership controls before dependency acquisition or R work.
2. Windows normally builds and installs the same checkout. It runs the entire installed test suite, standalone preservation and registered-binding gates, native C++ controls, parser failures/width/precedence, VariantAnnotation deferral, VCF write/index/query/sibling/error-cleanup controls, and four small scientific editions. The latter use the existing public R worker and complete independent19-field golds, source/audit/checksum/RDS checks and cross-format equality. Windows then runs all five legacy operations and compares complete decoded outputs/types with that run's Linux gold.
3. **Installed package evidence** starts after the producer terminates. Its receiver comes from the trusted workflow-definition commit. It fetches independent run/attempt/jobs/artifact service records and the exact original ZIPs, validates source/build/installed/runtime/library identities and all results, and records its own decision. Downloaded artifact tools and scripts are never executed or imported.

A producer's green aggregate is not independent terminal acceptance. Acceptance requires complete installed-suite coverage, the native/scientific gates, and direct five-operation Linux/Windows parity. Historical fallback assertions for Windows do not satisfy the last requirement. Old Windows or Linux cohorts do not qualify a newly built checkout.

## Source and dependencies

Git's actual commit/tree/blob identities and tracked file bytes are frozen before building. R builds a new source tarball from an exact work copy. Every included R/native/installed-test file must match this checkout; old, missing or additional code/test files reject. The generated archive, ordinary install command and actual loaded namespace/DLL are retained. Installed tests are copied byte-for-byte into owned scratch and are never rewritten by CI.

The public closure pins each source URL, version, compressed length/SHA256 and archive expansion length. If the exact CRAN source leaves src/contrib with HTTP404/410, acquisition may try only that identical file under the official same-host Archive/package directory. Both attempts/statuses are retained; a hash mismatch, other HTTP failure, newer version or changed Bioconductor release never supplies a replacement. Linux installs148 code packages plus the two GRCh37 reference packages actually used by its fixtures (150 total). Windows installs the same148 code versions plus all six references needed by its complete suite (154 total). Used reference source archives and complete installed extdata payloads must match across hosts. Native binaries need not match across operating systems. Historical seqminer is absent.

Package, dependency and base-library inventories must match before and after the gates. Namespace versions use R's normalized version spelling when compared with raw DESCRIPTION versions; physical paths remain exact. Actual R/compiler/runtime and Linux system-package identities are retained.

## Ordinary tests and development

ALS tests use distinct, byte-verified private copies of the shipped fixture. Their generated BGZF/index/query evidence stays outside the installed package. The tilde test uses a fresh Rscript child with both HOME and R_USER established before startup and keeps its original assertions.

The helper also supports a genuine pkgload development namespace with an already-compiled source DLL: it uses documented load_all(compile=FALSE, helpers=FALSE) and checks source/DLL identity. An installed CI run cannot fall back to this mode. Scratch cleanup is scoped to the calling test. An optional existing MSS_TEST_ARTIFACTS directory receives a verified evidence copy before scratch is removed, including failure evidence.

## Resource and failure evidence

Windows retains the accepted owned Job Object/no-breakaway boundary, hard7 GiB Job commitment limit,512 MiB global commitment reserve and original job-start-based300/305/316/317-minute compute/drain/seal/final deadlines. The Jobs API start time counts checkout/setup time too. Linux uses the existing pidfd/subreaper ownership/drain algorithm and reports a sampled aggregate7 GiB RSS cap; this is explicitly different from Windows' hard Job commitment cap.

Both lanes admit actual memory and free disk before acquisition and work, retain the source expansion/build bounds and keep8 GB free. The current conservative initial free-disk gates are16,996,512,195 B for Linux and66,420,158,684 B for Windows. Runner documentation is not an actual capacity observation. An unavailable source, changed digest, failed capacity gate, native control failure, timeout, skip, missing assertion or scientific difference fails qualification and retains available evidence. No custom runner or administrator setup is assumed.

The upload uses a source-sealed physical snapshot after confirmed owned-process closure. On missing closure/seal, only bounded named diagnostic prefixes are copied; a live tree is not recursively uploaded. Success and failure envelopes stay distinct. The terminal receiver preserves original service records/ZIPs and does not turn a partial or failed run into a pass.

## Running and receiving

Use the producer's manual workflow dispatch, or the configured push/PR paths. For a PR, retain its head SHA, tested merge SHA and merge parents separately. The receiver workflow has read-only contents/actions permissions and no personal or cloud credentials.

The evidence workflow can also be dispatched for an explicit completed producer run and attempt. Its source CLI is:

    python -I -S -B tools/ci/installed/service.py terminal --repository OWNER/REPO --run-id RUN --attempt ATTEMPT --out FRESH_SERVICE_DIRECTORY
    python -I -S -B tools/ci/installed/receiver.py --service FRESH_SERVICE_DIRECTORY --out FRESH_RECEIVER_DIRECTORY

The first command uses optional locally configured GITHUB_TOKEN only at GitHub's API. A signed artifact redirect is followed without that token. Signed URLs are not retained in logs. The receiver handles artifact data in Python and needs no R/compiler. Retain the receiver's trusted source revision with its receipt.

Initial publication still requires reviewing the full source integration and a genuine hosted run. The source proposal and metadata controls alone supply no new Windows qualification or performance measurement.
