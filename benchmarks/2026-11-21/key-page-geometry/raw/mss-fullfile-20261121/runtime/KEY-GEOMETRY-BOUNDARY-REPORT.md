# Typed key boundary source handoff

The sibling driver appends three format-3 lookup cases to `parity` after the existing block checks. The executed counter driver remains at its original path and pin. The sibling preserves the original matched mode, pressure mode, parity fixtures and assertion guard. The reviewed trial header patch remains the sole cache geometry change. [Boundary brief](KEY-GEOMETRY-BOUNDARY-BRIEF.md)

## Artifacts and pins

| Artifact | SHA-256 |
|---|---|
| [Executed counter driver](/private/tmp/mss-key-geometry-experiment-20261121/source/lookup_parity.cpp) | `505b0e8c72a517e7db4b48fb49126c0bfcbeeb6254dcd19c8090db57dd20fec9` |
| [Sibling boundary driver](/private/tmp/mss-key-geometry-experiment-20261121/source/lookup_parity-boundaries.cpp) | `b11f9d1897aacc32fc3899f546f27c6958a908a7df648f6bc1a5c24d37cf7888` |
| [Boundary-only review diff](KEY-GEOMETRY-BOUNDARY.diff) | `16a2c3d949b3c6013be9a6171c54bc1d3d50547adb657f891d687da4c65ddf6b` |
| [One-line trial header patch](/private/tmp/mss-key-geometry-experiment-20261121/source/trial-header.patch) | `cdfc635afd0c1b7f654301801c12d0f5bbab4393da830b9b0e259ebef81c96cf` |

The diff uses the original scratch filename as its logical path. The boundary source is staged separately under the sibling filename.

## Fixture coverage

Each fixture contains two sorted key rows, the existing tiny FASTA and its FAI, and one typed `A` record `A\t1\t1\tA\tC\tr\t1\n`. The driver writes exact little-endian key-row offsets with the existing `write64` helper and queries through ordinary `Reference::lookup`.

| Case | Zero-based key-file layout | Assertions and transcript |
|---|---|---|
| Final-byte newline | `I:` plus 16,379 `a` bytes and `\t0\n` place LF at 16,383. `I:r\t0\n` starts at 16,384. Offsets are 0 and 16,384. | Both aliases return one known candidate and consume one unit of work at cap 1. Both complete candidates print through the existing `fields` helper. |
| Malformed next-row offset | The same first row is followed by `I:r\t01\n` at 16,384. Offsets are 0 and 16,384. | Queries for valid key `I:r` with work 0 and 1 require the exact typed integer error and preserve the corresponding work value. The transcript prints both initial and remaining work. |
| Split UTF-8 identifier | `I:` plus 16,381 `a` bytes places `雪` octets `E9 9B AA` at 16,383–16,385. The first row ends with `\t0\n` and the second starts at 16,389. Offsets are 0 and 16,389. | The driver asserts the exact octets, runs two long-short query pairs, and checks and prints both complete candidates on each pair. |

Successful lookups check chromosome `1`, position and original position `1`, REF `A`, ALT `C`, ID `r`, the singleton identifier and original-position sets, and the empty original REF and allele fields of the typed `A` record. The malformed case requires `Typed reference has an invalid unsigned integer`. A two-row binary search visits the second row at its first midpoint. `visit_entry` parses the typed offset before invoking the visitor, while matching context work is consumed inside the later lookup visitor. The work-0 and work-1 checks characterize that ordering. [Typed reader and lookup](/private/tmp/mss-lazy-variantannotation-20261121/src/preservation_core.h:395)

## Source verification

- Checked the executed source and trial patch SHA-256 pins. Removing the new helper block and the single parity-only call from the sibling reconstructed the executed driver byte-for-byte. The read-only `git apply --check --whitespace=error-all` passed against the original scratch driver, and the diff reports 78 insertions and zero deletions.
- Independently evaluated all three source byte layouts, row sorting and concatenation, the UTF-8 octets, and both little-endian offset entries per fixture. Checked that the actual fixture code writes offset 0 and the computed first-row length.
- Reviewed the successful candidate assertions against the typed `A` decoder and the malformed case against `typed_uint`, lower-bound traversal and context charging. Checked the order of the alternating UTF-8 calls and the exact existing error text. The committed product checkout at `f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b` showed no tracked changes.

Baseline and trial builds, assertions, and byte-for-byte parity transcript comparison remain pending qualified guarded execution. Public API timing and promotion remain later gates.
