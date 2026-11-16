# Sealed count coordinator corrections — independent review

**Specification verdict: RED. Quality verdict: changes required.** The concrete
original failures are corrected, but an Important identity-alias/ordering case
remains within the original placement/seal-gate finding. No full-reference count
or production build is authorized. No fixes were applied by this reviewer.

Reviewed HEAD `2b80bbe72c01ff11edb74ba1fb2b3e4641f16bb4`, behavioral source
`0674d09930be48541d5e5306a0e72ebcb1164fd8`, based on preserved
`547c7d46004c38f6a1f545748027af911ef28542`. A new no-hardlinks detached clone at
`/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-fixes-stage11`
remains clean. Coordinator SHA256 is
`6d152cbd926f8949821a055f709023cecb8f20329b7f5873630b43c3b6fc069c`.
The complete 150,732-byte correction diff was saved/read before judgment; its
index portion was fully parsed and every explicit artifact verified. Original
requirements, the RED findings and the independently supplied root-context ruling
were included in the assessment.

## Remaining Important finding

**P2 — A source-plan hardlink to counted data is opened before both guard seals.**
[sealed_count.py:236](/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-fixes-stage11/sealed_count.py:236)
checks source-plan and producer-guard control locations lexically after rejecting
symlinks, but does not apply the hardlink identity refusal used for external
profile/admission controls. The per-build loop then fully processes the first
guard, source plan and receipt before checking the second guard.

`probe_alias_gate.py` constructs only a newly owned tiny fixture. It hardlinks
`GRCh37/variants.tsv` to the external `outside-source-plan.json`, supplies that
alias as the malformed profile's source-plan pin, and adjusts the synthetic first
guard command consistently. The second guard remains `running`. Calling
`validate_seals` hashes and parses the alias (two recorded `rb` opens of the same
device/inode as the counted TSV), then refuses with `Invalid/truncated JSON`.
Thus the counted file is content-opened before the other producer's guard is
accepted. Eventual JSON refusal is too late for the explicit no-input-open gate.
The checked control has link count 2, so metadata-only refusal is available.

This is a deliberately malformed local profile and synthetic guard context. It
does **not** imply that a real producer succeeded with a TSV manifest, that a
reviewed source plan had this identity, or that the program must sandbox malicious
trusted Python. It tests the same malformed-control/placement defense as the
original publication-contained-artifact regression, now through an identity alias.
No actual producer or active input was used by this probe.

Evidence: `probe_alias_gate.py`, `hardlinked-plan-gate.json`, and its owned fixture.
Reproduction uses the existing Python environment and a fresh owned output root;
do not overwrite the preserved evidence directory. The result is deterministic
and does not need process-metadata access or a native child.

A separate normal-control ordering probe records first guard → first source plan
→ first receipt → second guard (three bounded opens each), followed by refusal of
the running second guard. Both guard status/command checks can precede source-plan
and receipt content work using the already supplied expected roots and guard pins.
Artifact content hashes already occur after both full seal passes. Guard-first
ordering alone cannot protect a guard file which is itself a counted-data alias;
small external controls still need metadata identity protection before content
opens. A conservative single-link requirement for those controls, or explicit
device/inode comparison against the known counted files, is sufficient in scope;
no filesystem discovery framework is needed.

Counted data may legitimately have hardlinks. A separate positive links a counted
TSV externally, retains genuine distinct source-plan/guard files, completes the
worker and leaves publication filenames/file digests unchanged. This confirms that
the fix need not reject hardlinked scientific data. See
`guard-order-and-data-link.json` and `probe_guard_order_and_data_link.py`.

## Corrected behavior and positive evidence

The independent expected-root API has no profile-inferred fallback. Public and
internal entry points reject missing context, publication-contained case/profile/
admission paths and symlink parents before controls open or the first mkdir.
Parsed roots must equal the independent context before source seals are opened.
The original directory-write and profile-placement reproductions now refuse;
guarded actual-tiny success leaves both full publication trees unchanged.

Expected and observed 16MiB JSON limits are checked before initial and repeated
hash/parse reads. Exact-cap JSON remains valid, and bounded reads retain digest
and device/inode/size/time/mode stability checks. Necessary population, record/key
byte, maxima, ALT-disposition, rejected-header and FAI framing checks reject the
original impossible-response injection before planning success. Sixty independent
real-native responses cover empty input, invalid/overflowing positions, empty ALT,
unsupported/multiallelic/repeated ALT, Unicode and duplicate aliases, CRLF, empty
fields and missing final newline; all were accepted. Twelve independently generated
malformed neighbors refused. This validates malformed-response defenses, not a
claim that the pinned real native binary previously emitted contradictory values.

The C++ source has no correction diff. Its exact original 20 numeric fields match
the original diagnostic on both newly produced tiny publications. Their 3 emitted
contexts each include the historical alias expansion. Two separate tiny writer
oracles match records/keys/offset/rejected bytes and bound FAI bytes. No production
R/src/API/default, producer, guard or scientific predicate was changed. The capacity
formulas, 8GiB/160GiB limits, finite deadline/output caps, 170GiB planning reserve
and unconditional no-admission fields remain intact.

Fresh scoped verification:

- 35 selected correction/core tests passed; 67 unrelated cases were deselected.
- 60 independent real-native valid neighbors passed; 12 malformed neighbors refused.
- One actual tiny producer/writer/guard integration passed, with 20-field baseline
  parity on both builds, exactly 2 tiny writer oracles, and unchanged publications.
- The new source-plan-alias probe reproduced the remaining early-open defect.
- A legitimate hardlinked-data positive passed without source mutation.

The author's frozen 105-case log and selected receipts were hash-verified, not
rerun wholesale. Unchanged C++ sanitizers were not repeated. Known ordinary
counter/writer/original-counter binaries were copied after hash verification;
there was no native compilation, environment installation or lock change.

## Preservation and limits

All 272 author index entries verified, totaling 2,347,154 bytes. Author index SHA256:
`42ad13072e0a5ca1549d594ce9e474afb749009f0c67b7734c879600c51c869d`.
Both compressed large proofs roundtrip exactly and match their original raw
size/digest bindings. The copied original probe scripts are byte-identical; AST
inspection confirms exactly one occurrence of each of the two relocated path
literals in each script. The replay harness changes those literals in memory;
its executions are not represented as byte-identical executions at original paths.
Baseline source/helper snapshots equal the exact original commit.

The original RED report/index pins remain
`034b2b1e80089cb934f7112c26b50d001424ae6856c6bd47a45e72fe8cad8313` and
`7df58bbb6973abc31cb756b5aef9a1eada508cd0cd69c298f544c6e7e7e46ff2`.
All 555 entries of that original review index were reverified unchanged, including
its original source and tiny evidence. No old evidence path was reused. The
preliminary FAI correction and earlier writer-oracle totals retain their original
scope; this review adds only its 2 tiny writer oracles and zero production builds.

Both real reference preparations remain unsealed/out of scope. No active TSV,
SQLite, FASTA, archive or staging contents were opened, hashed, scanned or changed;
no active process was signaled. There was no new prefix/full count, full index,
network/provider/resource/credential/public/OAuth action, broad task-base discovery,
memory access or HOME/home/CODEX_HOME reassignment. This is no timing, peak-memory,
Windows or whole-file qualification. Every successful worker still reports
`full_reference_admitted=false`, zero count-process/production index builds and
`decision=no_admission`. A future GREEN correction review could permit only
consideration of a separately admitted full-count profile, never execute it.

Owned evidence and the explicit `evidence-index.json` are under
`/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-fixes-independent-stage11/`.
The index pins this report, complete diff, scripts, known clone sources/binaries
and owned tiny proofs. It excludes itself, symlinks and unrelated data.
