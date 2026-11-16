# Sealed count coordinator — scoped alias/order correction

**Decision: no admission; scoped independent re-review pending.** New no-hardlinks
clone/branch:
`/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-alias-fix-stage11`,
`fix/sealed-control-alias-stage11`. Behavioral source is frozen at
`08d93bf4b8d010b8ae73bf20663ec41c975d73c2`, based on preserved
`2b80bbe72c01ff11edb74ba1fb2b3e4641f16bb4`. Root alone may later admit a separately
pinned full-count profile; this task performs no full count or admission.

The remaining Important review finding was reproduced before changing source:
“A source-plan hardlink to counted data is opened before both guard seals.” The
normal-control trace also confirmed that the first plan/receipt was processed
before checking the second guard. The sealed review/report index pins are
`cda633dc61194719ca531a5dd2ed600b6e44b532f624868ff9579530f654d550` and
`046af9c7f33b708eae69d654a68c2578582306140f6459553e8a2597281e1fd9`.

## Focused change

Small external source-plan and producer-guard files must now have one filesystem
link during metadata-only preflight, before any control content hash/parse. This
closes the demonstrated counted-data identity alias, including a guard-as-alias
neighbor. It intentionally does not impose a single-link rule on scientific
TSV/FASTA/contig-alias inputs: legitimate counted-data hardlinks remain valid and
no huge scientific file copy is required.

Validation now authenticates **both** guards' successful status and exact expected
commands before reading either source plan, producer receipt or artifact content.
Publication/receipt path checks in that first phase are metadata only. The second
phase retains all existing source-plan/receipt semantics and output identities;
artifact hashes still follow complete seal validation. Counted-input hashing still
waits for both complete authenticated seals. The independently supplied expected
publication roots remain mandatory before any profile/admission/case content or
mkdir. Existing bounded JSON, identity/digest, capacity/model and resource gates
are unchanged; no C++/scientific/MSS API/default or guard implementation changed.

The malformed synthetic alias profile is not evidence that a real producer
succeeded with a TSV manifest. This is the same placement/identity refusal class,
not a malicious-Python sandbox requirement or a discovery framework.

## Observed RED→GREEN and focused verification

Both supplied probe scripts were copied byte-for-byte. An explicit harness changes
only their two hardcoded ROOT/CLONE path literals in memory via AST, using separate
new RED/GREEN destinations. Original review outputs were not reused. The copies'
SHA256s remain:

- `probe_alias_gate.py`: `be6a22dd6ff1efeab5b5f7f676271b3879a2d94302f1458786a074ab41ada7bf`
- `probe_guard_order_and_data_link.py`: `0556c3d49dde853e7c991cb4cabe048aa85ff0ecbc4d422b2ef62883cc4b8756`

Before the fix, the exact source-plan alias probe recorded 2 `rb` opens of the same
counted TSV identity, then invalid-JSON refusal while the second guard remained
running. Afterward it records **0 counted-identity opens**, refusing on small-control
hardlink metadata. The normal running-second-guard trace changes from first guard
→ first plan → first receipt → second guard to **first guard → second guard only**.
There are no plan/receipt content opens in that failure path. The supplied valid
hardlinked-data positive succeeds before/after and leaves publication files intact.
Historical replay receipt source paths point into this new clone before edits;
a byte-identical baseline source snapshot is pinned for auditing against 2b80bbe,
rather than treating the later corrected source at that path as the baseline.

The 12 new tests observed **8 failures/4 passing positives**, then **12 passed**.
They cover source-plan/guard aliases for both builds, running/wrong-command second
guards, both-success guard-first order, all three hardlinked counted input types,
and the conservative single-link rule on small controls alone.

The final appropriate subset plus one actual tiny producer/writer positive passed
**36 cases**: 12 new alias/order cases, 23 existing placement/context/size/stability
cases and 1 actual integration. Exact argv, inherited environment plus the single
`MSS_BASELINE_COUNTER` path override, output location and exit 0 are recorded in
`alias-evidence/verification-command.json`. The command uses the existing Python,
`-m pytest -vv --basetemp=alias-evidence/final-focused`, and the explicit listed
node selectors; it does not run the unchanged broad 105 suite.

Both newly owned actual tiny publications count all 3 emitted contexts each and
match all 20 original numeric fields from the pinned original diagnostic. Exactly
**2 tiny native writer oracles** independently match records/keys/offset/rejection
bytes and bound FAI bytes. Their full publication filenames/hashes remain unchanged
by counting. The actual owned guard succeeds naturally; no censor tests or process
signals occur in this pass. The three direct hardlinked-data positives retain
identical device/inode links and content, with no copy or count-process index.

Every success still reports `decision=no_admission`,
`full_reference_admitted=false` and `index_builds_executed=0`; these index fields
refer to count/production behavior, not the 2 tiny writer test oracles. No full,
prefix or current real reference is opened/counted. No unchanged C++ compilation,
sanitizer replay, performance/memory-peak/Windows qualification, environment/lock
installation, network/provider/resource/credential/public/OAuth action occurs.
Only exact known source/control paths and newly owned tiny fixtures are used;
there is no broad filesystem/evidence walk. No active TSV/SQLite/FASTA/archive/
staging content is read/hashed/modified, and no active process is signaled. No
.remember/tool memory is read/copied; HOME/home/CODEX_HOME remain inherited.

## Preservation and pins

Both old repositories remain clean at 547c7d4 and 2b80bbe. Their named report/index/
source pins and the sealed review pins were rechecked, as were all 13 preliminary
FAI artifacts. Old fixture trees, reports/indices and the original 39-file archive
are not copied or modified; this is not another whole-index evidence walk. The
historical 76 FAI writer-oracle correction and prior ordinary/sanitizer totals keep
their original scope. This pass adds only its 2 actual tiny writer oracles and 0
production/count-process index outputs.

| Artifact | SHA256 |
| --- | --- |
| Corrected coordinator | 094ca45bc80e5576c303ad392f5836f3297a5e48234f69ec81aad99c8789de50 |
| Alias/order tests | 53262ce8aecf9714dc2cf86fb3f829fe7bfbc7e80edc8d7f647c547c6e5332d5 |
| Replay harness | 193520629231b511809eea06bc2e3f107c356005a63e6da61e4fb4b6d60d7088 |
| Unchanged native source | fb529bff15666905d39c5726d2f49638a15c8ee61d90880da1b13079a954a10b |
| Unchanged ordinary counter | fa133e2c65bf9ed4ea0cb52791293b19b3a01fbead95f089ae13c0dff2b6c7fd |

`coordinator-alias-fix-evidence-index.json` is the compact explicit allowlist for
source, copied probes, observed RED/GREEN logs and traces, selected linked-data
identities, actual tiny producer/writer/parity evidence, preservation/command
receipts and this report. It excludes itself, symlinks, redundant unit trees,
large raw boundary-test fixtures, inherited historical indices already bound by
preservation hashes, unrelated/active data and memory. No broad fixture copy is
performed. An evidence-only commit may follow the behavioral freeze without source
changes. Scoped re-review and any later full-count admission remain outstanding.
