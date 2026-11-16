# Reviewed count coordinator — v1/v4 and v2 producer provenance

This is operator-only counting and logical planning. Exact reviewed source:
ca476ba84ebd0841d0409e06f08e69d7f912c2da, coordinator SHA256
963308cb023b1b729e53f891390db717f86dd2132853d1122063b013c3013e70.
The bounded independent review passes; it does not provide full-count readiness,
scientific equivalence, execution permission or an index/resource action.

Legacy profile v1 keeps seven global artifact roles and v4-only behavior. Profile
v2 has six common roles and REQUIRED per-build preparation_version/producer_source
pins (exact path,size_bytes,sha256 shape). For v2, only the two reviewed v4/v5 byte pairs
in dependencies.json are accepted. No global fallback, redundant global producer,
source import/version inference or v4 receipt relabeling is permitted. Authentic
mixed v4/v5 seals need not be re-prepared solely for operator bookkeeping.
Trusted completed producer context must authenticate historical source/command/
receipt association: a hash now cannot prove what a past process imported, and
trusted entrypoint argv need not literally name the module path.

Independent expected publication roots remain required outside profile JSON,
before profile/admission/case content or mkdir. Both guard statuses/exact commands
precede plan/receipt/source content; both complete producer seals precede source
hashes and counted input. Small controls, including v2 source paths, have metadata
placement/symlink/single-link protection. Scientific TSV/FASTA/alias hardlinks remain
allowed without copying. Sources are read-only; owned output cases stay outside
both publications. Known source size/JSON16MiB/numeric bounds and input/source
pre/post identity/digests remain mandatory.

Every success permanently says decision=no_admission, full_reference_admitted=false,
index_builds_executed=0. This NEW coordinator requires a NEW full-profile/source-
commit review and explicit admission if any future full count is considered.
NONE is provided or generated here. Old permissions/profiles are not rewritten.
The current root run is not inspected, altered, adopted or authorized by packaging.

The unchanged sampler covers the synchronous known worker/native-child lifetime;
it is not a new production parent-exit/orphan-aware supervisor proof. 8GiB/160GiB,
finite deadline/output caps, 170GiB reserve and both-order math are unchanged.
No new performance, peak-memory, Windows or whole-file qualification is claimed.

Runtime remains the exact coordinator/counter source plus unchanged guard_common.py
and process_memory.py. Native binaries and the 4aff... MSS header remain externally
pinned, not vendored. Frozen v4/v5 module snapshots live only in the public evidence
archive as .txt source proofs, never as extra runtime adapters. COORDINATOR.md is
exact reviewed v2 documentation. REVIEWED-README.md remains historical source-era
context; the new author/review reports and dependencies.json define this update.
The prior 94-file stage14 archive is preserved byte-for-byte, with its original
76/42/40/etc tiny-writer scopes and original 1bf 2.864s prefix-only timing attribution.
This v2 verification phase ran 0 new writers/compilers/full counts; packaging runs
no tests, native/producer/compiler/bench/data/provider operation or admission.
