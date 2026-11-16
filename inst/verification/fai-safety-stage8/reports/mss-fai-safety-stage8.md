# FAI numeric safety: public coverage and frozen qualification

**Implementation qualification GREEN; independent review, full installed-package
testing and actual native multi-OS CI remain required before adoption.** No
throughput benchmark was run or relabelled. No cloud operation, frozen-library
replacement or full-reference preparation was performed.

Frozen candidate: **`03d11e92b417add0183dea4fb77716f5502d5735`**, clean branch
`feat/fai-safety-validation` at
`/private/tmp/gwas-harvest-stage2/mss-fai-safety-stage8`.
Production fix remains the coordinator's
`6d6ec4187c69e113eda28fd2d68d986254a330c4`; this follow-up adds public/portable
tests, four native adverse cases, one native valid terminator control and contract
documentation. Delivered baseline is
`392d8266062fcc843c415cd0894914a79b147a33`. The initial plan/native RED and all
subsequent source states remain in Git and the original results directory.

## Scope and code review

FAI fields and consistent sequence-line termination follow the
[HTSlib faidx specification](https://www.htslib.org/doc/faidx.html).
The change checks numeric parsing, names, widths and addressed byte ranges. It
does not independently reconstruct or authenticate every FASTA byte after an
operator deliberately rewrites a trusted manifest and FAI together.

The checked unsigned decimal parser permits leading zeros and refuses empty
values, signs, whitespace, suffixes and uint64 overflow. FAI names must be unique,
nonempty and free of NUL/whitespace. A nonempty contig requires positive line
bases, no fewer bytes than bases, and a final addressed base inside the actual
FASTA file. Empty contigs preserve the native writer's zero-width representation
without entering division. An empty FAI is refused. FAI CRLF and an unterminated
last index row remain supported by the reader's existing line handling.

Critical arithmetic review found no additional in-scope production defect:

1. For `z=size-1`, the last addressed base is
   `offset + floor(z/bases)*bytes + (z mod bases)`.
2. Before multiplication, the code checks
   `rows <= available/bytes`. Nonempty-layout checks guarantee `bytes>=bases>0`.
   Therefore `rows*bytes` cannot overflow or exceed `available`; the subsequent
   subtraction and tail comparison are safe. `bytes>=bases` makes base-address
   order monotone, so bounding the final base also bounds earlier requested bases.
3. `sequence()` rejects nonpositive positions before subtracting one. It checks
   `zero<=size` and `n<=size-zero` before reserving memory or doing line division,
   avoiding the old unsigned addition wrap. Positive requests to empty contigs
   fail before division; zero-length requests keep their existing empty result.

There is no restriction to one- or two-byte sequence terminators. Numeric safety
does not imply that an intentionally re-pinned index's names/offsets describe the
intended biological sequence; that separate authentication claim is excluded.
Arbitrary-size malicious FAI text is not claimed to have a new constant memory
limit. Existing immutable-artifact checksums and operator resource guards remain.

## Meaningful RED evidence

The coordinator's original `red.log` records **17 of 20 invalid constructors
accepted**, three old refusals, and six valid controls. That log and initial source
were retained unchanged.

The new public R tests construct tiny format-2/3 references, then deliberately
re-pin each malformed FAI's MD5 and size in its fixture manifest. Consequently the
normal full checksum gate succeeds, and rejection must come from the native FAI
reader. They never query through a malformed accepted baseline index.

Corrected public RED against frozen baseline `392d826` has **46 assertion failures
and zero errors**: 21 of 24 malformed cases are accepted per format; two additional
old numeric exceptions per format have no FAI-specific diagnostic. All **156 valid
control assertions pass on the same old code**. This distinguishes intended safety
refusals from valid-input regression. `public-red-v2.log/.csv` is the qualifying
RED. Initial `public-red-v1` also caught a test-only incomplete column mapping;
that setup error was diagnosed and corrected, not treated as a package defect.

The four additional malformed cases are trailing numeric whitespace, hexadecimal
suffix, an empty-contig offset beyond the file and an embedded NUL name. The
native NUL fixture uses explicit string concatenation; an initial test-only
manual string-length mistake was corrected before the final native/sanitizer run.
Historical logs are retained and superseded by final pinned source/evidence.

## Verification actually executed

| Gate | Final result |
|---|---|
| Focused public FAI tests | **2 blocks / 300 assertions**, no failures/errors |
| Re-pinned malformed public FAI cases | 24 cases × 2 formats; every rejection has native `FASTA index` diagnostic |
| Valid public layout controls | LF, CRLF, multiline/single-line no-final-newline, empty contig, UTF-8 name; both formats |
| Decimal representation controls | Leading-zero fields and CRLF FAI text; exact scientific row parity |
| Actual SourceCpp portable entrypoint | **82 blocks / 1,727 assertions**, no failures/errors/warnings |
| Actual registered `.Call` entrypoint | **42 blocks / 814 assertions**, no failures/errors/warnings |
| Complete native script | Core, block I/O, projection, typed-format structural/occurrence/caps and FAI fixtures pass |
| Focused native FAI | **24 invalid refusals / 7 valid controls**, range boundaries pass |
| AddressSanitizer + UndefinedBehaviorSanitizer | Same native FAI/range fixture passes; no sanitizer findings |
| Mathematical SNP corpus, each format | **16,128 rows**, zero violations |
| Mathematical haplotype/contract corpus, each format | **79,904 rows**, zero violations |
| Explicit scientific/resource corpus, each format | **297 rows**: 156 resolved, 132 nullable, 9 refusals |
| Complete native output parity to reviewed scientific 7 | **96,032 rows per format byte-exact**, all 12 native fields |
| Existing `open_reference` Rd and diff checks | Passed; no R API/Rd change |

Public positives compare exact source, all 19 normalized columns, independently
expected coordinates/ref/alt/effect/EAF and unchanged normalizer/version fields.
Leading-zero re-pinned indexes are also opened through the explicit relocation API.
The native valid controls cover LF, CRLF, no EOF delimiter, leading zeros, Unicode,
writer-generated empty contigs and a consistent vertical-tab sequence terminator.
Maximum `size_t` requests at positions 1/2 and a maximum signed position are
refused before allocation; ordinary indexed sequence remains exact.

The finite mathematical oracles use independent expectations, then additionally
compare complete outputs with the previous reviewed scientific-7 format2/3
artifacts. The SNP corpus retains 4,608 resolved and 11,520 empty uncertainty rows
per format. The haplotype corpus independently validates 2,104 resolved payloads
and 17,736 paired ALT–ALT explanations, with no violations. No new resolution,
allele orientation or sign improvement is claimed.

The new R file is loaded explicitly by **both**
`tests/standalone/preservation.R` and `tests/standalone/registered.R`, including
the actual portable workflow's SourceCpp entrypoint and the registered bridge.
The latter builds `RcppExports.cpp` with the native DLL using only Rcpp/testthat
and the ordinary compiler toolchain. These are local macOS runs; actual Windows
and other-OS CI for this candidate is a coordinator gate, not an inferred pass.

## Guarded execution and identities

The scientific corpora and both R entrypoints ran under the unchanged **8GiB
memory and 160GiB free-space guard**: 23.07s qualification duration,
613,384,192-byte sampled peak guard memory, 469,473,067,008 bytes minimum free.
The focused sanitizer run independently passed the same guards: 1.68s,
240,451,584-byte sampled peak, 469,358,690,304 bytes minimum free. These are test
execution receipts, not throughput benchmarks. Guards monitor only the launched
process tree; no active preparation was signalled or changed.

| Artifact | SHA-256 |
|---|---|
| Production `src/preservation_core.h` | `4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc` |
| New public FAI tests | `8188ca8c72aeb11e198570e1a3ea5acc70012d38e4fdfb8c06f31939e8dd834d` |
| Expanded native FAI tests | `f528a15563cf65463088ff87f6c1b8621fed7c75bb2216b48c3600dfe3ff473c` |
| SourceCpp runner | `b8e4ffce12163c17e98d8706112c79dfbb5cff297d7199cefca1f3172310e5cf` |
| Registered runner | `0708ff9797654ae1ca0fa1ca85f7c380ded7284a81b80fc33e65701c352c19d8` |

R wrappers, generated registration files and public APIs are byte-identical to
delivered `392d826`. The complete suffix from `inline Normalized normalize(` to
the end of the core is byte-identical to that baseline, SHA-256
`59ab9d077a3d2fe3e7809bcdd566505ea7a30073120bd7c172c7338560020d67`.
Writer formats, scientific algorithm/version identities and valid output remain
unchanged. Older measurements keep their older explicit source/DLL pins; none
were relabelled with this new core hash.

Evidence: `results/mss-fai-safety-stage8/source-pins-final.json`,
`candidate.patch`, `qualification-summary.json`, per-format oracle outputs,
`qualification-guard.json`, `sanitizer-guard.json`, `native-final-v2.log`,
`public-red-v2.log/.csv`, `public-green-v1.log/.csv`, `Rd-check.log` and
`evidence-index.json`. Initial root native RED/GREEN and intermediate author logs
remain separate. Full installed suite, independent review and actual native CI
have not been run by this implementer follow-up.
