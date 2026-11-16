# Independent final MSS scale packaging review — stage 17

Spec verdict: **PASS**. Packaging quality verdict: **PASS**. The timing-comparison finding is **addressed**; no new blocking or nonblocking finding was identified in the fix.

Reviewed Harvest `d495faaa7f317502e46922bed00325cb5e65ba89` and MSS `3670586851168920c8c183e7d34479b76199f92d`, against package commits `d48bc80d75cd7ed4925220b603f93d53601629dc` / `a3eef80df652fd4211c416b804bda82889496090` and bases `d997af5dd15fdaadd0ba573795e326e4a798abff` / `df8870104868a0e3d946a889f6aff6beb6aac468`.

## Timing correction

Both preserved reports contain the same scope error: source-result SHA256 `1fced5ffd819520771e5bf728f064be9f87fa24908306f8bc052a5a045a7e498` and packaging REPORT SHA256 `730b401ee6a798a32cc5c0f93245c584c315582dc53abf4afeb8e35a2bbbd831` mix final processing 224.500 seconds with legacy external wall 142.250 seconds. Erratum SHA256 `34ff9e82b5b9ce933c1246a9ded1406bc24c457e507cb00a0b63ecb43bdfce84` explicitly supersedes this comparison in **both** reports. It gives external wall **229.370 versus 142.250**, processing **224.500 versus 137.107**, and peak R RSS **1.857620239 versus 16.943 GiB**. Neither timing comparison shows a final-over-legacy speedup. Different APIs/scientific/output work, single shared-host observations, uncontrolled OS caches and repeated-SNP scope are explicit.

I read and independently hashed the existing legacy metadata `harvest/docs/verification/stage2/mss-bench.md`: 14,143 bytes, SHA256 `5d66dc6297703ed7ab9ffb109b54c7c44500c9bbfafacad77a81ab8baa68a06a`. Its provenance identifies adapted legacy `b4a2d21`; the original matched table supports all three legacy values. Its “cold” label expressly excludes flushed OS caches. Both current status pages lead with the authoritative erratum and correct wall/processing comparisons. Both erratum copies equal the external addendum. The supplied erratum diffs agree with the actual committed additions and exact status-page insertion.

## Immutable package checks

Each archive directory has exactly **47 regular files / 139,099 bytes**: 46 explicitly indexed payloads plus the index. Each indexed size and SHA256 matches; the replicas are identical. Index SHA256 is `53ead16a1621846eeea7e6187c25767441bafdd295cfef2e3fe1aa2990d38e16`. Both 27,866-byte gzip tar copies hash to `9293227d361c1477830d75d547747902d47fe8dbee4b50d151d4a6b380cf3339`; their 184,320-byte raw tar hashes to `b0bb156f365570ebb28b0dcc035d9bad73760a0eb0bc9c9b2cea61a28d02c9ff`. All 47 regular members match the explicit file set and directory bytes, with no duplicate or extra member. All eight gzip logs decompress to the bound exact original small logs, with matching raw sizes/digests and zero gzip mtime. Both release JSON copies and the external release agree.

Commit/tree comparisons independently establish 50 documentation-only additions in each original package commit, no modification of any pre-existing tracked path, and exactly two subsequent erratum changes: the new erratum and status-page insertion. All 51 resulting added files are Git tracked and their working bytes match their committed blobs. All 47 archive blobs, tar and release are unchanged by the fix. The pinned old94/counter40 trees and tar/release blobs retain their recorded identities. Existing source, R/src, production and main-plan tracked paths are unchanged at the reviewed commits. No unrelated integration change was introduced.

## Receipt and scope checks

I compared 32 explicitly selected historical metadata files to their original retained bytes, including the original report, all selected run summaries/correctness/metrics/specs, 10M/40M guards and fixed completion/operator/pin receipts, root acceptances/authorities and separate counter-review metadata. All match. Both author evidence indices and every indexed file's size/digest match. The 10M/40M original failed-stage paths are absent; exact root receipts require operator exit zero plus accepted fixed completion and no failed stage. No target was rerun.

All three final metric receipts name measured production `03d11e92b417add0183dea4fb77716f5502d5735` and DLL SHA256 `6f23605e06ac0b35fd40268b08190fd55e9a09e2fccded2f8cd7b629a3cd49eb`. Every row has a passing all-19-field receipt, all rows examined, two exact source copies and zero unaccounted deltas.

| Rows | Processing s | External R wall s | Peak R RSS GiB |
| ---: | ---: | ---: | ---: |
| 1,000,000 | 24.444 | 29.010 | 1.5177612305 (1.518 rounded) |
| 10,000,000 | 224.500 | 229.370 | 1.8576202393 |
| 40,000,000 | 811.171 | 815.070 | 2.2537231445 |

The package distinguishes processing, external R wall, oracle verification, whole operator and aggregate logical storage. It describes a repeated 437,458-record real chr22 biallelic non-palindromic SNP panel with synthetic exact beta/EAF, without a 40M-unique or complete-study claim. It preserves the historical legacy40M manual-censorship limitation. External driver/control/source/library/reference anchors are explicitly external; no new source/DLL/scientific-payload verification is implied. Historical authority is archival. Sampled synchronous guards do not establish hard instantaneous ceilings or general orphan-aware supervision. Full reference preparation/count/index, full corpus and new-source admission remain outstanding.

This review used only bounded explicit metadata/log/archive reads, digests, gzip decoding and Git metadata/object identity comparison. It executed no shell/Git subprocess, target tests/imports, R/native/PG/producer/compiler/writer, provider/network action, process inspection, timing rerun or scientific payload walk. It read no .remember/tool memory, private configuration, credentials or authentication logs. Only this review's own report and evidence outputs were written. Permanent disposition remains **decision=no_admission; full_reference_admitted=false; index_builds_executed=0**.
