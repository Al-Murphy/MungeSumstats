# Native reference count-only refinement

**Candidate, pending independent review.** No production input or native index is
admitted by this diagnostic. Source is frozen at
`1bf1a08685488798496e9d7978b1792dbd290a9d`, clean branch
`feat/reference-count-only` in `results/reference-count-only-stage9`.

This implements the exact dry-run option requested by the earlier typed3 capacity
report. It reads a plain FASTA and prepared TSV, reuses the unchanged MSS normalizer
and emits numeric byte/count JSON. It writes no native records, sort runs, keys,
offsets, FAI or index directory. There is no full-input launcher, public route,
build permission or MSS R/API/default change. Production header SHA256 remains
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`.

The retained payload count is exact for successful native emission. Every record
has a distinct physical byte offset. Since a complete sort row includes that
pointer, equality across records is impossible; per-record unique keys determine
final key bytes/counts without sorting. Raw sorter counts preserve duplicated
original/canonical coordinate keys and repeated primary/alias IDs. Unsupported
contexts and invalid positions follow the original writer's disposition. Repeated
historical context rows remain separate physical occurrences. Checked uint64
arithmetic refuses overflow. No source row is dropped to reduce a bound.

## Verification

- Initial stub:24 expected differential failures, retained in baseline-red.log.
- Current29 cases pass in0.36s, comparing actual native writer records/rejections/
  keys/offset file sizes for formats2/3. Cases include indels, multi/unsupported/
  empty ALTs, duplicate/high-fanout aliases, repeated contexts, missing IDs,
  Unicode, invalid/signed/overflow positions, CRLF and no final newline.
- The duplicate-alias fixture independently checks19 raw versus10 final key
  rows and exact raw byte multiplicity using the real sorted key file.
- A malformed diagnostic version suffix first reproduced acceptance, then was
  fixed to require exactly2 or3. The RED receipt is preserved. This is not an MSS
  product parser change.
- Standalone overflow controls prove integer addition and ostream counting refuse
  overflow without wrapping the prior count.
- Separate ASan/UBSan builds pass the unchanged29 assertions and overflow fixture.
  A fresh strict run enables halt-on-error and explicitly checks empty stderr for
  all56 successful native invocations plus overflow control. The earlier runner's
  successful-subprocess stderr inspection gap is documented and its evidence kept.
  This is macOS/Apple compiler evidence; native Windows diagnostic parity is not
  claimed. Original production source remains unchanged.

## One frozen prefix calibration

The existing31,408,384-byte /1M-site prefix and3,339,671,278-byte FASTA were pinned
and hashed before/after, along with the production header, executable and harness.
The source SHA is
`3136a1b6bb0dd7e6bba90410bacf47ac2e10eaab0d07e190a38a1be8ed066372`.
No existing index was rebuilt and no active source preparation was read.

| Quantity | Exact observed / prior pinned expected |
| --- | ---: |
| Source rows / S records |1,000,000 /1,000,000|
| Successful A records / raw ALT tokens |1,242,071 /1,242,077|
| Rejected ALT occurrences |6|
| R, record bytes |63,883,154|
| U, raw sort bytes |100,705,522|
| M, raw key rows |4,726,213|
| K, final key bytes |75,230,231|
| N, final key rows |3,484,369|
| Offset bytes |27,874,952|
| J, rejected bytes |167|

All values exactly match the earlier actual typed3 build/recount. Maximum observed
source line/record/key-row lengths are8,663/8,650/24bytes. This prefix contains no
historical alias tokens and is incomplete/nonrepresentative. The synthetic cases
exercise aliases separately; neither predicts full history-expanded inputs.

The guarded **count process and descendants** complete in2.864s, sampled peak RSS
6,389,760B /physical5,292,416B, minimum free457,121,832,960B. Its guard limits are
8GiB process-tree/160GiB free/16MiB owned receipt case/600s. Source/FASTA files are
existing external pinned inputs already occupying the measured filesystem, with
their3,371,079,662bytes explicitly declared. The measured phase excludes the
coordinator's streaming pre/post hashing and is not an MSS harmonization speed
measurement or a claim that all operator memory fits6.4MB. Case receipt size at
sampling peaks3,606B; no index scratch is produced.

With an explicit25,906-byte FAI +19-byte marker +1MiB manifest allowance and no
contig-alias file for this prefix, the retained logical upper size is199,471,389B;
the conservative phase envelope is268,768,527B. The excess over the actual prefix
retained result is the explicit manifest allowance, not a changed record model.

The envelope is the maximum of retained size, data plus two coexisting raw-sort
copies, and data plus last raw run plus final keys/offsets. Allocated filesystem
overhead, directory entries, simultaneous work and170GiB planning reserve remain
additional. Both sealed input closures/build orders need independent sizing and
fresh free-space measurements. These numbers authorize no production build.

## Remaining gates

Independent source/proof review is required. A future full-input coordinator must
authenticate complete prepared receipts, declare its source/FASTA/core/executable
pins, bounds and small-file allowance, guard all phases, verify unchanged inputs
and retain failures. Both complete GRCh37/38 preparations are still active and
were neither scanned, changed, signaled, copied nor cleaned by this work. No full
count, full native index, cloud offload/resource, billing or reference default
change occurred. The user's delivered MSS/Harvest feature branches are unchanged
by this isolated candidate.
