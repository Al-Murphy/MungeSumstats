# Timing comparison erratum and packaging addendum

This erratum is authoritative for the legacy 10M timing comparison. It supersedes only that comparison in both frozen reports:

- Original source-result report `reports/mss-final-source-scale-stage17.md`, SHA256 `1fced5ffd819520771e5bf728f064be9f87fa24908306f8bc052a5a045a7e498`, including its exact copy inside the 47-file archive.
- Packaging `stage17-final-mss-scale/REPORT.md` and its exact external copy `reports/final-scale-package-stage17.md`, SHA256 `730b401ee6a798a32cc5c0f93245c584c315582dc53abf4afeb8e35a2bbbd831`.

The frozen wording compares final 10M processing time 224.5 seconds with legacy 142.25 seconds. **142.25 seconds is legacy external wall time, not processing time. Those timing scopes must not be compared.** The corresponding processing value is 137.107 seconds.

| 10M implementation | Processing seconds | External R wall seconds | Peak R RSS GiB |
| --- | ---: | ---: | ---: |
| Historical adapted legacy `b4a2d21`, original cold matched workload | 137.107 | 142.250 | 16.943 |
| Final production `03d11e92b417add0183dea4fb77716f5502d5735` | 224.500 | 229.370 | 1.857620239 |

The valid wall comparison is **229.370 versus 142.250 seconds**; the valid processing comparison is **224.500 versus 137.107 seconds**. Neither comparison shows a speedup of the final implementation over legacy. The memory comparison is **1.857620239 versus 16.943 GiB**, scoped to this repeated-SNP workload. APIs, scientific/output work and implementation versions differ; these are single shared-host observations with uncontrolled OS cache, not matched-equivalence performance ratios or cold-disk medians. The older metadata uses “cold” for initial application/reference state and explicitly says OS caches were not flushed. Its legacy source is recorded as `b4a2d21`, separately from its older new/final preservation implementations; the measured new production source here remains exact `03d11e92b417add0183dea4fb77716f5502d5735` and DLL `6f23605e06ac0b35fd40268b08190fd55e9a09e2fccded2f8cd7b629a3cd49eb`.

The authoritative retained legacy metadata is `/private/tmp/gwas-harvest-stage2/harvest/docs/verification/stage2/mss-bench.md`: **14143 bytes, SHA256 `5d66dc6297703ed7ab9ffb109b54c7c44500c9bbfafacad77a81ab8baa68a06a`**. Its opening provenance identifies the older adapted legacy source, and its “Original matched measurements” table records `legacy-10M-cold` as external wall 142.25 seconds, processing 137.107 seconds and peak RSS 16.943 GiB (lines 60–78 in the reviewed file). This is a documentation-only metadata anchor; no historical benchmark or scientific payload was opened or rehashed.

The original reports, all 47 archive files, their index and compressed/raw tar remain byte-identical. Index SHA256 remains `53ead16a1621846eeea7e6187c25767441bafdd295cfef2e3fe1aa2990d38e16`; compressed tar SHA256 remains `9293227d361c1477830d75d547747902d47fe8dbee4b50d151d4a6b380cf3339`. Source numeric receipts, individual run times and all 1M/10M/40M correctness results are unchanged; no test or benchmark was repeated. This separate addendum corrects the reporting label and comparison, not measured data or execution history.

Read this addendum with the [frozen package report](stage17-final-mss-scale/REPORT.md) and [exact release metadata](stage17-final-mss-scale-release.json). It grants no new benchmark, full reference count/index, profile/admission, deployment or provider authority. Root owns the combined main-plan link; the separate status page leads to this erratum.
