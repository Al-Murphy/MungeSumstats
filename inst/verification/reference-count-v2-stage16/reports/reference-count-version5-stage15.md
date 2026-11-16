# Explicit per-build producer provenance for count profile schema2

Implemented in the new isolated clone
`/private/tmp/gwas-harvest-stage2/results/reference-count-version5-stage15`.
Base: `306aea51ea45eb8ae781b339029d2b1b93b5388f`.
Head: `ca476ba84ebd0841d0409e06f08e69d7f912c2da`.
Code/tests/docs are committed and the tracked/untracked working tree is clean;
owned evidence is excluded locally and explicitly indexed. The original stage11
clone, earlier author/review evidence, delivered runtime/archive and active
producer/source publications were not modified. Root independent review remains
required before packaging/adoption or any separate admission.

Schema1 retains its seven global artifact roles and v4-only behavior. Schema2
requires six common artifact roles and a mandatory exact `producer_source` pin per
build. No global fallback, redundant global producer, version inference/import,
or relabeling of v4 receipts is permitted. The two reviewed compatibility pairs
are fixed, not calculated from convenient replacement bytes:

| Version | Frozen Git provenance | Bytes | SHA256 |
| --- | --- | ---: | --- |
| biodb-reference-stream-4 | 84199d3b37273ab56c838aac385ffcdea40197bc | 38909 | c6327b05066a8a00508c724cd90ac44a0df2ccaa86da43831b43d888dd2b9fb5 |
| biodb-reference-stream-5 | 6b30709918a8bf73f5ff51dd1be5483e96ddc39d; same reviewed module in release979fd325 | 43423 | 7a7668c91235ddf8037563ae3141b04e8fd75604c1f24348198a858fb5186147 |

The only frozen module acquisitions executed these exact argv, with stdout
written into the clone's `fixture_sources/producer_v4.py` and `producer_v5.py`:

```
["git", "-C", "/private/tmp/gwas-harvest-stage2/biodb-history-empty-stage14", "show", "84199d3b37273ab56c838aac385ffcdea40197bc:src/biodb/gwas_reference_stream.py"]
["git", "-C", "/private/tmp/gwas-harvest-stage2/biodb-history-empty-stage14", "show", "6b30709918a8bf73f5ff51dd1be5483e96ddc39d:src/biodb/gwas_reference_stream.py"]
```

No active producer source file was opened. The ordinary existing native binary
was copied from the original stage11 clone, without compiling: 88816 bytes,
SHA256 `fa133e2c65bf9ed4ea0cb52791293b19b3a01fbead95f089ae13c0dff2b6c7fd`.
Its source remains unchanged (8595 bytes,
`fb529bff15666905d39c5726d2f49638a15c8ee61d90880da1b13079a954a10b`), as does
the reviewed header (49600 bytes,
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`).
`implementation-pins.json` also records the unchanged guard/harness pins. These
are identity/provenance checks, not new compiler/native qualification.

Both source paths join metadata-only placement/symlink/single-link checks.
Declared source version/size/SHA metadata must match its reviewed pair. Both guard
status/exact-command checks precede all plan/receipt/source content. Both complete
seals precede producer-source hashes and counted input content/hashes/native calls.
Known actual source size is checked before every source hash, including repeated
snapshots; all source pins appear in pre/post control snapshots and per-build
results. Mutation during the worker refuses success. Trusted completed producer
context still must authenticate historical imports; current hashes alone do not.
An exact trusted argv may use an entrypoint that imports the module without
naming its path as an argument.

## TDD and executed scoped verification

The binding brief and read-only preflight/required-tests matrix were read first.
TDD and verification skills were applied; the approved ruling supplied design.
Tests were written before implementation. Initial RED:

- `pytest -q test_producer_versions.py --basetemp=version5-evidence/red-fixtures`:
  **9 failed, 33 passed**, with the missing schema2 gate causing the three
  positives, five guard-order cases and worker-mutation case to fail. Initial
  negative cases already refused under unsupported-schema behavior; the GREEN
  run exercises their intended provenance/placement paths.
- `pytest -q test_actual_versions.py --basetemp=version5-evidence/actual-red-ps-fixtures`:
  **2 failed, 1 passed**: actual schema1/v4 counted successfully, actual v5 and
  mixed profiles refused at the old schema gate.

The first attempted `python3 -m pytest -q test_producer_versions.py
--basetemp=version5-evidence/red-fixtures` could not start tests because default
Python3.15 lacks pytest; subsequent commands used `/opt/homebrew/bin/pytest`
(Python3.14). The first actual RED command used
`--basetemp=version5-evidence/actual-red-fixtures` and had **3 failures** from
sandbox refusal of `ps`, before any producer child spawn. The bounded local
process-metadata exception was then used for actual guard tests only; no approval
review rejected an action.

Exact subsequent pytest argv, relative to the new clone, and observed results:

```
pytest -q test_producer_versions.py --basetemp=version5-evidence/green-fixtures
# 42 passed
pytest -q test_actual_versions.py --basetemp=version5-evidence/actual-green-fixtures
# 3 passed
pytest -q test_producer_versions.py test_coordinator.py test_coordinator_corrections.py test_control_alias_gate.py --basetemp=version5-evidence/regression-fixtures
# 163 passed
pytest -q test_actual_versions.py test_coordinator_integration.py -k 'frozen_preparers or censor' --basetemp=version5-evidence/final-actual-fixtures
# 5 passed, 1 deselected
pytest -q test_producer_versions.py test_coordinator.py test_coordinator_corrections.py test_control_alias_gate.py --basetemp=version5-evidence/final-regression-fixtures
# 163 passed
```

Final verification is **168 scoped cases passed** (163 focused regressions plus
5 actual/resource cases), including 52 newly authored cases. `git diff --check`
passed. The unrelated broad application/native writer/sanitizer suites were not
run, per the brief. The deselected old producer/writer integration is explicitly
outside this stage; its historical result is preserved.

Coverage includes wrong/missing/unknown/wrong-type/swapped receipt/expected
versions; schema1/v5 rejection; missing/wrong/swapped source pins; missing common
roles; global-only/redundant producer authority; wrong declared and actual source
sizes (including an oversized sparse owned source); same-size content change;
source mutation in worker; inside-publication/symlink/hardlink/data-alias controls;
both guards first and complete seals first; legitimate shared source references;
trusted argv without literal module path; schema2 repeated source/JSON bounds;
old alias/order/data-hardlink/read-only publication cases; exact16MiB JSON,
numeric/bool/addressability/native count invariants, emitted-row equality,
full-admission refusal, resource deadline/output censor, and no-admission behavior.
No native C++/writer/arithmetic/FAI logic changed.

Actual positives use complete newly owned tiny FASTA, VCF (full multi-ALT, symbolic
ALT, unsupported star and indel context), assembly, same-release merged/withdrawn
history and rights receipts, under the unchanged synchronous process-tree guard:

| Final actual pair | Schema | GRCh37 version/rows | GRCh38 version/rows |
| --- | ---: | --- | --- |
| legacy | 1 | v4 / 3 | v4 / 3 |
| empty-history plus valid one-to-many | 2 | v5 / 4 | v5 / 4 |
| mixed authentic source pins | 2 | v4 / 3 | v5 / 4 |

V5 receipts record one empty merge record and complete supplied history; the
nonempty route fans out to both live rsIDs. Each count equals that build's trusted
expected emitted rows. Both final source/publication snapshots remain identical;
no `fasta.fai` is created in a publication.

Across the actual RED-with-ps, initial GREEN and final actual runs, separately
counted activity is **18 tiny producer executions**, **7 successful guarded count
workers / 14 native counts of actual tiny publications**, plus **2 synthetic
censored workers**. The sandbox-blocked attempt spawned no producer. Additional
ordinary tiny counter subprocesses occur in focused unit regressions; those are
not represented as a grand total by the actual-publication counters. This stage
executes **0 new tiny writer oracles, 0 compilers, 0 full-reference counts and 0
index builds**. Historical writer counts (including the earlier76 FAI correction
and stage11's2 writer oracles) retain their original scope and are not reset.
`executed-argv.json` contains 41 explicit producer/guarded-worker/native/censor
argv records; native argv are the exact frozen coordinator construction plus
successful worker results. Logs and selected206 tiny control/output files are
explicitly pinned, without recursive evidence discovery.

## Scope and evidence

Permanent results remain `decision=no_admission`,
`full_reference_admitted=false`, `index_builds_executed=0`. The 8GiB/160GiB limits,
finite output/deadline, 170GiB reserve and both-order planning math are unchanged.
No real/full profile or admission was created or executed, no full preparation,
restart, active payload scan/hash, benchmark, provider, native writer/compiler,
subagent, push/deployment or user-repository write occurred. No `.remember`
contents were read, copied, deleted or indexed; inherited HOME/home/CODEX_HOME
were not repurposed. No broad filesystem/evidence walk occurred. Sparse oversized
source/JSON files belong only to new negative fixtures and were refused before
content reads. This does not establish scientific equivalence of v4/v5, historical
code-execution attestation, new orphan supervision, peak-memory/speed guarantees,
Windows qualification, or full-count/index admission.

The explicit index is
`/private/tmp/gwas-harvest-stage2/results/reference-count-version5-stage15/version5-evidence-index.json`.
It pins this report, the requirements/preflight receipts, explicit code/log/control
lists and selected newly owned tiny evidence. Full metadata/argv/counter receipts
are in `version5-evidence/`. Independent root review and any later admission are
separate work.
