# Independent review: Harvest worker relocation and MSS portability fixtures

**Verdict: GREEN for both scoped candidates. No blocking findings.** This review
does not change production source, installed libraries, cloud resources or full
reference preparations.

| Candidate | Frozen review commit | Scope |
|---|---|---|
| Harvest | `95c55d663eb05e0820668173bcf1890f7d606a19` over `b6c1f07` | Select the SHA-256-pinned runtime FASTA and pass an explicit MSS override only when needed |
| MSS | `50912b79d93ac7ee14127cb3a83fb4c16a5060d5` | CI registered bridge plus portable relocation fixture corrections |

Fresh `git archive` exports were used in
`results/worker-relocation-independent-review-stage7/{harvest,mss}`. All **918
Harvest files** and **768 MSS files** remained byte-exact after the runs. Author
checkouts remained clean. While the review ran, the coordinator advanced the MSS
checkout to archival-only descendant `74c85cdb3a5959d1460103212e8db336c01769bc`;
its `R`, `src`, `tests` and `.github` are unchanged from the frozen review commit.
The independent archive remained pinned to `50912b7` throughout.

## Actual local verification

| Run | Result |
|---|---|
| Installed MSS relocation-capable worker profile | **13 tests passed**, no skips/errors/failures, 134.47 seconds |
| Installed earlier scientific7 worker profile | **13 tests passed**, no skips/errors/failures, 132.26 seconds |
| Fresh corrected MSS SourceCpp suite | **80 blocks / 1,427 assertions**, exit 0 |
| Fresh corrected MSS registered `.Call` suite | **40 blocks / 514 assertions**, exit 0 |
| Deterministically denied timestamp-write fixture branch | **8 blocks / 192 assertions**, no failures/errors/warnings |

Each worker profile ran the author's eight actual integration tests plus five
independent configuration/closure tests. The older-library profile asserts refusal
of relocated references; it does not claim that older MSS supports relocation.
Both profiles used actual installed R/C++ workers and tiny local references. DLL
hashes were rechecked before and after:

- Relocation-capable: `b63e4e43a167ad5b4807bc6a36758e820eb7b0c70ae3e77124e6ab909d9e4841`.
- Earlier scientific7: `c1c0a8b29479c33d8b4cfa4646e46e39c4da436686f4ca2da12a402251003d50`.

The independent worker runs used process-tree guards with a 4 GiB memory ceiling,
160 GiB free-disk floor, 2 GiB owned-evidence ceiling and ten-minute deadline.
Initial sandboxed attempts stopped before starting tests because `ps` was denied;
the same bounded local harness then ran with approved process-monitoring access.
The guard was retained, not bypassed. No cloud/network operation was performed.

## Harvest worker assessment

The Python model already enforces a unique `fasta` label and unique file paths.
The new configuration obtains the physical FASTA from that validated member after
`ReferenceBundle.verify()` checks its SHA-256. R reads the unchanged RDS manifest,
normalizes its recorded path and the supplied runtime path, and either preserves
the old one-argument API or requires the optional `fasta` formal. MSS independently
checks original FASTA size/MD5 and reference artifacts. Python then compares the
complete live R handle closure with the declared bundle before recording a ready
worker receipt.

The supplied and independent tests establish:

- Same-path use works with both installed libraries, preserves one native handle
  across jobs, and adds no relocation fields to normal edition metadata.
- Transferred references remove the original paths, preserve manifest bytes and
  content identity, and yield exact source artifacts, source TSV, reference
  checksums and all 19 normalized fields. Two jobs reuse one verified native handle.
- A different, explicitly pinned byte-identical FASTA is used by upgraded MSS even
  when the recorded old FASTA still exists. The old library refuses that transfer
  rather than falling back to the available recorded file.
- A stale configured old path refuses before a scientific job. A configured
  byte-identical FASTA that differs from the selected bundle member also refuses:
  upgraded MSS opens that file, then Python rejects the live closure mismatch;
  older MSS rejects the unsupported explicit override during startup.
- A byte-identical `records` pin outside the actual index directory passes byte
  verification but fails live closure validation. This verifies that matching
  content does not silently substitute a different operational file set.
- Wrong FASTA bytes with stale SHA-256 refuse in Python. Deliberately repinning
  those wrong bytes still refuses against MSS's original manifest MD5. Neither
  refusal starts a scientific job.
- Every independent startup refusal leaves no worker receipt, `job-*` request or
  `measure-*` request. A zero-row version health probe remains the existing startup
  check and uses no scientific source.
- Physical paths remain in operational configuration/receipts. The scientific
  content identity remains unchanged on a byte-identical transfer.

No RDS rewrite, guessed path, weakened integrity check, new scientific algorithm or
public activation is introduced. The earlier ID/allele uncertainty, mutation,
version, crash, existing-output and expansion controls were rerun in both profiles.

## MSS portability correction assessment

All production `R/` and `src/` files are byte-identical to the previously reviewed
`8334f1c9052ae266dc9df80a35bca07868af1423` relocation candidate. Changes after that
candidate are workflow coverage, test fixtures and archived evidence.

The timestamp fixture now checks the actual outcome of `Sys.setFileTime`. If the
write succeeds, the stat must change and subsequent reuse must refuse. If Windows
denies the write while a native FASTA stream is open, the stat must stay unchanged
and the same input must still produce the exact baseline result. The fixture then
opens a fresh binding and appends exactly one binary byte to unheld `source.tsv`;
it asserts the one-byte size change and requires the complete binding stat guard
to refuse. This preserves mutation coverage on every OS without treating a denied
mutation as a changed file. An independent run forcing timestamp writes to return
`FALSE` without changing bytes exercises that branch successfully.

Binary append avoids the Windows text newline conversion that invalidated the
earlier `+1 byte` assertion. `gc()` before offline directory transfers releases
unreachable native streams; it does not close a reachable handle or permit a live
reference to be modified. The final fixture uses normal copies/moves with UTF-8
paths, not privileged symlinks. The workflow now runs both SourceCpp and the actual
generated registration bridge on Windows, macOS and Linux with Rcpp/testthat and
ordinary compiler tooling; Bioconductor, WSL and external bcftools are unnecessary
for those CI steps.

The original failed run `43844881386` and diagnostic run `43847267264` were inspected.
The diagnostic log shows `set_time=FALSE` and unchanged before/after timestamps for
both reference formats. Later retained failures show the text newline byte-count
and open-file directory-move fixture issues. Their failure records are preserved;
production code was not changed to mask them.

## Independently checked coordinator gates

Saved native run **43848950725**, head `50912b79d93ac7ee14127cb3a83fb4c16a5060d5`,
has three completed successful jobs. I parsed the decoded logs using explicit
SourceCpp/registered phase boundaries: **each of Windows, Ubuntu and macOS ran
80 blocks / 1,427 assertions plus 40 blocks / 514 assertions**. The Unix standalone
C++ gold step also succeeded; its Windows omission is the existing intentional
workflow condition. The native run is
https://github.com/Al-Murphy/MungeSumstats/actions/runs/43848950725.

The coordinator-executed final installed-package gate was independently checked
against its CSV, summary, guard and **234 matching source pins**:
**255 blocks / 2,270 passing assertions, zero failures/errors**, two known warnings
and two missing-reference skips. Warnings are the existing a0/a1 and extra-rsID-info
tests; skips are the existing dbSNP-backed validation tests. These are not new
relocation failures. This review verified those retained results rather than
rerunning the memory-heavy installed suite. Earlier 2,262-assertion results remain
separate and unchanged.

## Evidence and reproduction

Evidence root: `results/worker-relocation-independent-review-stage7/`.
`final-summary.json`, profile `guard.json`/`junit.xml`/`run.log`/`dll-pin.json`,
`worker_probes.py`, `run-worker-review.py`, local MSS logs,
`denied-timestamp.R`/log/summary, `native-ci-verification.json`,
`full-installed-evidence-verification.json`, `initial-source-pins.json`,
`final-source-integrity.json`, and `author-evidence-verification.json` contain the
exact harnesses, commands, counts and hashes. All nine entries in the supplied
worker author evidence index were checked against their recorded sizes/hashes.
`evidence-index.json` pins this report and the independent evidence files.

From the evidence root, reproduce worker profiles with the prepared local Python
environment:

```sh
/private/tmp/gwas-harvest-stage2/harvest/.venv/bin/python run-worker-review.py supported
/private/tmp/gwas-harvest-stage2/harvest/.venv/bin/python run-worker-review.py unsupported
```

The harness refuses existing guard/log receipts; select a fresh output directory
before rerunning. The local MSS runners use `MSS_TEST_RLIB=/private/tmp/harvest-Rlib`
and the fresh `mss` archive. No inherited FAI semantic-hardening expansion or
cloud operation was included in this review.
