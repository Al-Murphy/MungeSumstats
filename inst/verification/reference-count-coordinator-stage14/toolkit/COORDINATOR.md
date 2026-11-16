# Operator sealed count coordinator

`sealed_count.py` is a local operator diagnostic. It never builds an MSS index.
The only public CLI runs the complete counting worker, including input hashing,
FASTA validation, both native counts and post-hashing, under the unchanged reviewed
process-tree guard. No HTTP service or JSON submission authority is provided.

A local operator must independently review and supply a profile, its expected
SHA256, both expected publication roots OUTSIDE profile JSON, and a new exclusive
output path. Missing independent root context fails closed; there is no inferred-
root fallback. All paths must be absolute, without symlinks or parent traversal.
Expected source/receipt hashes come from
the reviewed acquisition/preparation context. This tool does not construct a real
profile, discover a replacement digest, trust a producer's self-declaration alone,
or turn a newly observed digest into an expected one.

```
python sealed_count.py --profile /absolute/reviewed-profile.json \
  --expected-profile-sha256 EXPECTED_REVIEWED_PROFILE_SHA256 \
  --expected-grch37-publication /absolute/sealed-GRCh37 \
  --expected-grch38-publication /absolute/sealed-GRCh38 \
  --case /absolute/new-count-case
```

Before any control content open or first mkdir, metadata-only checks require the
profile, optional admission and output case to be outside both independently
expected publications. Symlinks/parents and publication directory aliases refuse;
hardlinked external controls refuse. Parsed final roots must then exactly equal
the independent context before seal/control processing. The internal worker and
test entry points require the same explicit context, and the launcher binds and
passes it into the guarded child; Python callers supply `expected_publications`.

The profile has schema_version1, explicit context (`tiny_fixture` or
`full_reference`), format_version3, and exactly GRCh37 then GRCh38 build entries.
See the newly authored fixture constructor in `test_coordinator.py` for a tiny
schema example; its synthetic guard/source artifacts are deliberately never full
producer authority. `test_coordinator_integration.py` additionally exercises real
tiny producer publications under the unchanged guard. Neither fixture is a real
full-source profile.

Each build pins its reviewed source-plan JSON, published `receipts.json`, completed
producer guard JSON and exact producer command, final directory, build/release,
preparation version and expected emitted rows. Both guards must have success,
exit_code0, reason null and natural_oom false. Their commands must equal the
trusted expected commands and contain the reviewed plan and final directory.
Every manifest field must match its receipt, including sources, reference policy,
scope, assembly and header assertions. Complete same-release merged and withdrawn
history is mandatory. Final directories, controls and three counted inputs must
be regular, nonsymlink paths. Pending publication directories refuse. Both guard
status/command checks finish before either source plan, producer
receipt or artifact is content-opened. Small external source-plan/producer-guard
controls require one filesystem link at metadata preflight, preventing counted-
data hardlink aliases before any hash/parse. This rule does not apply to counted
TSV/FASTA/alias data; legitimate scientific hardlinks remain supported without
copying them. Both complete seals are validated before opening any prepared
variants TSV, derived FASTA or alias TSV.

Full scope additionally requires `coordinator_source_commit` and a separately
root-reviewed, digest-pinned admission supplied with `--full-admission` and
`--expected-admission-sha256`. The admission must identify action
`sealed-full-count-only`, independent_review `GREEN`, the same source_commit and
the exact profile_sha256. The operator's local authenticated review is the trust
boundary; these JSON words are not signatures or external authentication.
No real full admission/profile was generated or used in this implementation.

Artifact roles are `coordinator`, `counter_source`, `native`, `header`, `guard`,
`harness` (the guard's process_memory.py) and `producer`. Each pin has absolute
path, size_bytes and lowercase SHA256. Their complete pre/post identity includes
device, inode, size, mtime, ctime and mode. All counted inputs bind expected size
and digest to the authenticated producer receipt. Input hashes and both counts
are sequential, without sampling, truncation, skipped rows or fanout caps. All
native numeric fields/maxima are retained and checked; source rows must equal
expected emitted_rows. Nonzero native status, successful stderr, missing/extra
fields, invalid JSON, booleans, inconsistent counts or changed pins refuse.
Expected AND observed JSON size must be at most16MiB before every hash/read,
including repeated external/control snapshots and parse reads. Exact-cap JSON is
accepted; stable file identity/digests still bind the complete parse. Read loops
never exceed the size observed before opening, and concurrent change refuses.

Native-response validation enforces necessary typed record/key byte minima,
zero/nonzero populations, maxima leaving room for every other row, ALT disposition
(and its source-byte bound), raw/final key/dedup relationships, raw key fanout,
rejected-header/status framing and positive FAI framing. An empty ALT string
produces zero native ALT tokens; invalid-position ALT tokens are not silently
classified as successful/rejected alleles. Empty, all-rejected, indel, multi-ALT,
Unicode and duplicate-alias output remains valid. These are malformed-response
defenses; the frozen pinned counter/producer remains the source of exact counts.

The profile must retain the 8GiB complete-tree memory ceiling and 160GiB hard
free floor, and supply a finite positive deadline and finite owned-output cap.
The unchanged guard supports macOS/Linux; it refuses when it cannot measure the
owned tree. This coordinator is not a Windows guard implementation or Windows
qualification. The native MSS/producer portability claim remains separate.
Censoring and measurement failure are failures, not successful counts or OOM proof.
The public CLI has no option for launching an unmonitored count worker.

For each build, C adds the native conservative FAI bound, an independently streamed
alias-table byte bound for R's quote-all/doubled-quote UTF8 output with conservative
CRLF framing, the19-byte format marker and the profile's explicit bounded RDS
manifest allowance/basis. The manifest bound is an operator planning assumption;
no RDS manifest is built or separately qualified here. Existing source bytes are
already occupied before actual free space is measured on the case filesystem.

```
I = T + R + K + 8*N + J + C
P = max(I, R + 2*U + J + C, R + U + K + 8*N + J + C)
A_then_B = max(P_A, I_A + P_B, I_A + I_B)
B_then_A = max(P_B, I_B + P_A, I_A + I_B)
```

All arithmetic uses exact nonnegative integers, checked uint64 sums and native
signed addressability. Allocation/directory/filesystem allowance, concurrent
reservations and170GiB planning reserve remain explicit separate fields. Equality
at the reserve boundary counts as a logical fit. Every result still says
`decision=no_admission`, `full_reference_admitted=false`,
`index_builds_executed=0`. These index fields describe count-process/production
behavior; verification intentionally creates separate tiny native writer oracles.
An actual build needs a separate native memory/max-line/sort-buffer profile,
independent review, fresh free space and post-first-build remeasurement.

Cases contain immutable stdout/stderr, worker pins/counts and the guard receipt;
`summary.json` appears only after guard success and post-exit resource/pin checks.
A failed/censored case must be preserved and a new path chosen for any retry.
Tests use internal Python helpers outside the CLI only on newly owned tiny fixtures.

Focused correction verification command (process-metadata permission is needed
for guard tests; the unchanged C++ sanitizer suite is not repeated):

```
/path/to/existing/python -m pytest test_coordinator.py \
  test_coordinator_corrections.py test_coordinator_integration.py -vv \
  --basetemp=/absolute/new-tiny-test-evidence
```

Sanitizer replay is `verify_coordinator_sanitizers.py`, after C++17 builds of
count_reference, writer_fixture and arithmetic_fixture into
`coordinator-evidence/sanitized` with address+undefined sanitizers and
`-fno-omit-frame-pointer`. It preserves all82 calls' raw stdout/stderr and refuses
any unexpected successful stderr. Full test receipts are indexed separately.

The original547c7d4 implementation/report/index and all preliminary FAI evidence
remain preserved. The independent stage11 RED findings are corrected only in this
new branch and still require fresh independent correction review. The byte-
preserved reviewer probes and explicitly logged path-only AST relocation harness
reproduce the original failures exclusively in newly owned tiny fixtures. Do not
rerun a probe over an existing evidence path. No full source is authorized here.

The scoped alias/order correction is a further isolated branch based on2b80bbe.
Its supplied review probes are preserved byte-for-byte with explicit path-only
relocation. Only the relevant placement/identity/order, size/stability checks and
one actual tiny producer positive are rerun; unchanged C++/sanitizers, broad105-
case verification and whole-reference work remain outside this correction.
