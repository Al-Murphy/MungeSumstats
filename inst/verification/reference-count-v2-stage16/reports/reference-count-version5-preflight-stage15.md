# Reference-count v5 compatibility — read-only preflight

Recommendation: retain immutable legacy **profile schema 1/v4** behavior and add a
small explicit **profile schema 2** for required per-build producer-source pins.
A profile-schema change is local operator bookkeeping: it does **not** require
re-preparing a successfully sealed v4 GRCh37 input, rewriting its receipt/plan, or
claiming it has v5 scientific semantics. Implementation and tests require a new
root ruling/dispatch; none are performed here. Full-count admission remains false.

## Confirmed blocker and identity gap

The frozen 306aea51 `sealed_count.py` accepts only schema_version1, demands the exact
7 global artifact roles including `producer`, and at its producer-version check
requires receipt version == per-build expected version ==
`biodb-reference-stream-4`. Thus genuine v5 output refuses even though v5 still
publishes `reference_format=mss-preservation-2` with the same six TSV columns,
output identities, emitted_rows and complete_supplied_history status. The v5
producer adds empty-history diagnostics/counts and semantic handling; those do not
require a native counter or byte-model change. Native format_version3, producer
preparation version4/5 and operator profile schema 1/2 are separate namespaces.

The global `artifacts.producer` is one code-file path/size/SHA pin. It is hashed
with other artifacts after both complete seal passes and in pre/post control
snapshots, but there is no build→producer-source binding. Guard commands are exact
per-build expected lists and must contain the plan and final directory; the
validator does not require the producer module path in the command. A single
file cannot honestly identify both distinct v4 and v5 module bodies for a mixed
pair. Merely widening the literal version allowlist or pointing the global pin
at v5 does not supply the missing v4 source association.

This is not a cryptographic proof of historical code execution. The trusted,
completed producer context must authenticate each guard/receipt/plan/source
association; hashing a source file now cannot establish that a past job imported
it. Do not require a literal module path in every production argv: a legitimate
entry script can import the module without its path being an argument. Continue
to compare the exact reviewed expected command and completed root provenance.

## Minimal bounded proposal

For v1 leave the existing v4-only roles/validation/format semantics intact and
keep frozen old source/evidence. Do not infer new mixed provenance from its global
pin. For v2 require the six common artifact roles (all current roles except
`producer`), and **each build** must supply `preparation_version` plus one
`producer_source` pin with existing path/size_bytes/sha256 structure. No global
fallback or simultaneous global/per-build authority is allowed. Both builds may
reference the same immutable source file when genuinely homogeneous.

For this rollout use the two reviewed pairs below, not versions inferred by
importing producer code or scanning strings. A two-entry compatibility check is
enough; no registry/service/framework is needed. Exact expected hashes come from
reviewed source/completed-job context, never freshly calculated replacements.
Code release 979fd325 and implementation 6b307099 have the same v5 module bytes;
checkout release HEAD must not be mistaken for a new preparation semantic version.

| Receipt preparation_version | Reviewed code provenance | Source bytes | SHA256 |
| --- | --- | ---: | --- |
| biodb-reference-stream-4 |84199d3 |38,909 |c6327b05066a8a00508c724cd90ac44a0df2ccaa86da43831b43d888dd2b9fb5 |
| biodb-reference-stream-5 |6b307099; documentation release 979fd325 |43,423 |7a7668c91235ddf8037563ae3141b04e8fd75604c1f24348198a858fb5186147 |

Validate declared version/source-pair metadata and pin structure without opening
scientific input. Put both per-build code paths into the existing metadata-only
outside-publication/alias preflight, and apply the small-control identity rule to
these new small source controls (shared references are fine; no huge data copy).
Known source size should be checked before hashing a wrong large path. Preserve
both guard status/command authentication **before** plan/receipt/producer-code
content. After both complete producer seals, pin actual source bytes with the
existing stable identity/digest logic. Include **both** source pins in every
before/after control snapshot, and record their versions/source bindings in each
build result. Source mutation must prevent the final success receipt.

The existing trusted-profile/admission hash and exact coordinator-source commit
binding still apply: a changed coordinator requires fresh independent review and
separate root admission, not reuse of old permission. No actual profile is created.

A same-schema 1 extension is technically possible, but only by defining two strict,
unambiguous authority shapes: the unchanged legacy all-v4 global form, and a fully
per-build form with no silent fallback or redundant global producer. That is a
real schema/meaning change despite retaining the number, and saves no meaningful
runtime work. A tiny v2 envelope is clearer than defaulting/overriding a single
global pin or retaining an unused global legacy artifact for homogeneous v5.
Neither approach justifies a 19-hour v4 re-preparation solely for bookkeeping.
The decision to accept an authentic v4 scientific preparation remains root's;
compatibility does not relabel v4 as v5 or assert equal biological coverage.

## Exact future verification needed

`required-tests.json` gives the scoped matrix, not executed tests. Required actual
positives are a frozen v4/v4 tiny legacy run, a real v5 tiny producer case including
empty history, and a real mixed 37-v4/38-v5 pair with separate immutable modules.
The current helper has one global PRODUCER and learns expected preparation_version
from the observed receipt; parameterize it by build and trusted expected source/
version constants instead. Preserve full ALT/indel/source context and compare each
counter source_rows with that build's authenticated emitted_rows; cross-build row
counts need not be equal. Existing pinned native binaries suffice if separately
authorized; no C++/writer/compiler change follows from this metadata compatibility.

Negatives must cover wrong/missing/unknown receipt/expected versions, missing/wrong/
swapped source pins, wrong actual source bytes/size, mutable sources, global-only
mixed authority and redundant global producer in v2. Assert **no counted input
content/hash/native call** before refusal and no successful summary. Also retain
both-guard-first order, aliased guard/source-plan refusal, new producer-control
identity/placement, no source-tree writes and allowed scientific data hardlinks.
The existing 16MiB JSON cap (including rehashes), numeric/bool/addressability checks,
emitted-row equality, complete owned synchronous worker guard, 8GiB/160GiB caps,
finite deadline/output cap, 170GiB reserve and both-order math remain unchanged.
No new orphan-supervision, Windows, peak-memory or speed claim is implied.

## Scope and receipts

Only narrow frozen validator/control/schema and actual tiny integration source
excerpts, the v5 receipt-control excerpt and the sealed public source-pin receipt
were read. Evidence records selected line ranges and their excerpt digests; the
full source 094ca45b/7a7668 pins are supplied reviewed identities, not a fresh full-
file scan. Delivered Harvest 9ffe equality is root-supplied, not independently
queried here. No active producer checkout/source path, full data, publication
receipt, archive/staging/SQLite/statistical payload or process status is inspected.
No code change, test/probe/app/native writer/compiler/preparation/count/full-profile/
provider execution, restart, signal, broad discovery, memory read or original
repository/evidence mutation occurs. Only this new report and receipt directory are retained as outputs; HOME/home/CODEX_HOME remain inherited. Root supplies the current running
v4/approved-v5 context; this report makes no new live-state or restart claim.

All permanent outcomes remain `decision=no_admission`,
`full_reference_admitted=false`, `index_builds_executed=0`. The preflight proposes
a bounded implementation; it neither authorizes that implementation nor any real
count. Explicit receipts/index are under `/private/tmp/gwas-harvest-stage2/results/reference-count-version5-preflight-stage15`.
