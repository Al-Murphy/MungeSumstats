# Verified runtime FASTA binding in the Harvest worker

Author candidate only; independent review is required before adoption.

Frozen source: `/private/tmp/gwas-harvest-stage2/harvest-worker-relocation-stage7`,
branch `feat/worker-reference-relocation`, implementation `95c55d6` over reviewed
Harvest `b6c1f07`. Plan `ff0415e`. Production changes are confined to
`mss_worker.py` and `mss_worker.R`; added integration tests are opt-in through an
explicit installed-library capability profile. No provider operation, RDS rewrite,
scientific algorithm change or public activation was performed.

Python sends the already SHA256-verified unique `fasta` member as operational TSV
configuration. R reads the unchanged manifest and compares normalized recorded and
runtime paths. Matching paths preserve the old one-argument open and older MSS
compatibility. A different path requires the reviewed optional `fasta` argument;
unsupported libraries fail during startup before scientific jobs. MSS still
independently verifies original FASTA size/MD5 and the entire index closure.
Reference membership from the live R handle must match the pinned Python bundle.
Physical paths remain outside the scientific reference content identity.

Actual installed-MSS regressions use tiny CPU fixtures and real R/C++ workers:

- Baseline RED: one transfer test failed with reference checksum/size mismatch
  after the original directory/FASTA were moved, proving the old worker opened
  the stale manifest path. Two cases were deselected, not run.
- Relocation library, production R source `8334f1c` and installed DLL
  `b63e4e43a167ad5b4807bc6a36758e820eb7b0c70ae3e77124e6ab909d9e4841`:
  **8 passed, 65.43seconds**, no skips.
- Frozen earlier scientific7 DLL
  `c1c0a8b29479c33d8b4cfa4646e46e39c4da436686f4ca2da12a402251003d50`:
  **8 passed, 62.55seconds**, no skips. This profile expects and asserts the
  unsupported-relocation startup refusal; it is not a relocation success claim.

The positive transfer removes original paths and uses spaces/Unicode in the new
location. Two scientific jobs reuse one native handle; all19 normalized fields,
source TSV/original bytes and reference checksums are byte-identical to baseline.
The RDS manifest bytes and scientific content identity are unchanged. Wrong FASTA
bytes are refused first by Python when SHA256 is stale; even with an explicitly
repinned wrong runtime SHA256, MSS refuses against its original manifest. Existing
worker mutation, version, crash, output and expansion controls run alongside.

Both commands use the unchanged process-tree guard with8GiB RSS/physical accounting,
160GiB free-disk floor,4GiB owned-evidence cap and600second deadline. Source/DLL
hashes were checked before/after. Guard/log/pin receipts and the exact harness are
under `results/worker-reference-relocation-stage7/`. No provider calls, reference
downloads, full production indexes or full-study harmonization were attempted.
Pilot orchestration still explicitly supports macOS/Linux; native MSS Windows
portability is a separate CI gate.
