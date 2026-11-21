# Matched namespace build producer review

**Specification compliance — changes required. Code quality — changes required.** Two Important findings affect resource enforcement. Neither requires a change to the frozen guard or the benchmark consumer.

| Artifact | Reviewed SHA-256 |
| --- | --- |
| Producer | `41da5092b5aafc65c7d95cbd684e53f4606e131b5f0a818293548080d502f4df` |
| Addition diff | `202c30487c5b91de5786db632423b3e64bb63cecb3ef6d994bb8afad97678fb2` |
| Implementer report | `7682b76dfc2fe8dc5730a3ee61fbdd5f55a97a9a883a0c2963eef1b6be8b8c1b` |

## Findings

1. **Important — a resource stop can trigger an unfunded write.** On every failure, line 299 passes `None` to the result writer. Lines 69–75 consequently open a new file and write up to 1 MiB without testing available space, owned usage or root identity. For example, if a guard phase reports that the owned directory has already crossed 512 MiB, this path adds the result to the same directory. If the filesystem is at its continuing floor, allocating even one block can cross it. The frozen guard explicitly refuses evidence it cannot fund. Initial admission rejection also reaches this writer, potentially with no room above the floor. Apply the existing funding and directory checks to the failure result when available. If evidence cannot be funded, exit nonzero, report the reason to stderr and preserve the existing evidence. A bounded size alone does not enforce the remaining budget. The same principle applies to the initial rejection record when free space is at or below the floor. [Failure dispatch](namespace-builds.py#L299) · [Writer](namespace-builds.py#L66) · [Admission](namespace-builds.py#L251) · [Frozen guard space accounting](/private/tmp/mss-performance-20261119/tools/darwin-full-stage185/guard.py#L673)
2. **Important — finalization can report completion after the shared deadline or a handled stop.** The final stop check precedes a full directory walk and provenance write. After that write, the operator leaves the signal context, marks the run completed and performs another full directory walk and result write. The guarded writer checks disk capacity only. If the deadline passes during either walk, the function can return zero with a completion record. A stop signal received while provenance is written can likewise be recorded by the handler and never checked. Keep success finalization inside the signal scope. Recheck the frozen guard after potentially lengthy final checks and before publishing success, and ensure a late stop makes the operator exit nonzero with an honest outcome. [Final stop check and provenance](namespace-builds.py#L285) · [Completion and exit](namespace-builds.py#L297) · [Writer](namespace-builds.py#L72) · [Frozen guard stop handling](/private/tmp/mss-performance-20261119/tools/darwin-full-stage185/guard.py#L540)

**Critical findings — none. Minor findings — none.**

## Consumer and execution scope

The successful-phase source uses the authenticated binding and its frozen guard, one admission time, the two requested commit IDs, direct Git archive stdout and separate private libraries. Both installs use the same explicit environment, SDK, thread controls and `en_US.UTF-8` locale. The runtime capture uses that same locale and checks the full captured row set against the original matrix. The producer retains actual install and runtime transcripts, terminal records, toolchain captures, archive and installed-file identities, native command tokens and staged-load markers in the shapes consumed by the binding. The top-level estimate and floor match the brief. I found no additional consumer-schema mismatch. [Role contract](namespace-bench.py#L123) · [Toolchain and runtime contract](namespace-bench.py#L191)

I checked the source, addition diff and implementer report pins. I reviewed the reported metadata controls without rerunning them. Native executables, package builds, payload hashing and benchmark measurements were outside this source review. Actual build success and independently observed terminals remain pending.

## Round 1 re-review

**Specification compliance — pass. Code quality — pass.** Both original Important findings are addressed. Open findings in the reviewed scope are 0 Critical, 0 Important and 0 Minor.

1. All operator JSON writes now use the frozen guard's current disk sample, funding arithmetic and validated directory handles. The guard binds the new root's device and inode before admission. A rejected admission can retain a funded rejection receipt. Failure reporting checks funding even after a stop, writes a diagnostic to stderr when a receipt cannot be funded and exits nonzero. [Writer](namespace-builds.py#L66) · [Root and admission](namespace-builds.py#L265) · [Failure reporting](namespace-builds.py#L310)
2. Both success controls are finalized inside the shared signal scope. The writer checks for a stop after the funding walk and after closing the new file. Its interruption cleanup operates relative to the validated parent, verifies the opened control's device and inode and unlinks only the matching name. A stopped completion attempt takes the nonzero failure path. A valid provenance file can precede that failure, so the documented handoff requires the completed result and zero exit. [Publication checks](namespace-builds.py#L70) · [Cleanup](namespace-builds.py#L81) · [Final result](namespace-builds.py#L303)

I reviewed the fix diff and the reported metadata controls. The producer, controls and fix diff matched the supplied SHA-256 values. Producer `615535da8374c81d3d7f23f0300c313d26d06ebea638e73239492243f851a65f`, controls `b757ed4a868217c3476a4b8c824671cb2e35116d7942e63f9d2193ef6e606b99`, fix diff `a232f116c4d60881deb6c595bc56b47e7aa9deca48a332ac22673563dcaff454`. The implementer report SHA-256 observed on disk was `df2aa847ba40843e6522666b877566a027630bac1b6772a231f4c4da6337f4f7`. Native execution remains outside this source re-review.
