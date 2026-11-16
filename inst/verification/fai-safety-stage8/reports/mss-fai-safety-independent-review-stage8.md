# Independent review: native FAI numeric and range safety

**Verdict: GREEN for the documented scientific contract of frozen `03d11e92b417add0183dea4fb77716f5502d5735`.**

Production core SHA-256:
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`.

No blocking finding. This verdict covers the documented native numeric/layout and
address-safety change and preservation of valid scientific behavior. It does not
claim authentication of an operator-rewritten FAI against all FASTA biological
content. Actual native multi-OS CI remains a coordinator adoption gate. No new
throughput benchmark was run or relabelled by this review.

Reviewed frozen candidate `03d11e92b417add0183dea4fb77716f5502d5735` against delivered
baseline `392d8266062fcc843c415cd0894914a79b147a33`, using fresh independent archive
exports in `results/mss-fai-safety-independent-review-stage8/`. Production change
`6d6ec4187c69e113eda28fd2d68d986254a330c4` is retained in candidate history. The
candidate adds tests/documentation after that production fix without modifying it.

## Independent execution

| Run | Result |
|---|---|
| Actual SourceCpp portable entrypoint | **82 blocks / 1,727 assertions**, exit 0 |
| Actual generated registered `.Call` entrypoint | **42 blocks / 814 assertions**, exit 0 |
| Complete native script | Core, bounded block I/O, projection, typed reference and FAI fixtures passed |
| Independent public R checks | **3 blocks / 86 assertions**, zero failures/errors/warnings, exit 0 |
| Independent big-integer FAI oracle | **11,271 cases**: 741 accepted / 10,530 refused; zero constructor mismatches |
| Exact queries through oracle-safe geometries | **25,935** expected-byte/refusal comparisons, zero errors |
| Physically consistent FASTA compatibility controls | **160 controls / 9,600 exact range comparisons**, zero mismatches |
| ASan + UBSan, both independent native programs | Both exit 0, zero sanitizer findings |

The native programs were independently written and compiled from the frozen core,
not copied from the implementation tests. SourceCpp and registration runners also
compiled the archive source freshly; no author-built DLL supplied these results.
The review ran locally on macOS arm64 with R 4.5.3 and local Rcpp/testthat.

The numeric corpus uses an independent Python arbitrary-precision last-address
formula rather than reproducing the implementation's division-before-multiplication
guard. It covers a systematic small-field Cartesian product, uint64 extrema,
3,000 deterministic random layouts, full-token numeric corruption, long leading
zeros, duplicate/empty/whitespace/NUL names, CRLF and unterminated FAI rows,
empty-contig EOF offsets, and maximum unused line strides. Query expectations use
Python's unbounded arithmetic and exact byte addresses, including maximum `size_t`
requests and maximum signed positions. Constructor cases intentionally test the
declared safety geometry; acceptance does not certify their biological mapping.

Physically consistent controls independently construct UTF-8-named sequences of
length 0, 1, 13 and 65,540 with line widths 1, 2, 7, 32,767 and 65,536. They exercise
LF, CRLF, vertical-tab and consistent two-byte alternative terminators, with and
without a final sequence/index delimiter. Requests cross the 65,536-byte page
boundary and include empty, negative, end, over-end and enormous lengths. Exact
expected DNA substrings are computed with widened arithmetic. This checks valid
compatibility separately from the geometry-only corpus.

The allowance for consistent non-LF terminators follows the
[HTSlib faidx specification](https://www.htslib.org/doc/faidx.html), which defines
LINEWIDTH as physical bytes per sequence line and allows consistent alternate line
termination. No extra one-/two-byte terminator restriction was introduced.

## Code assessment

The complete-token decimal parser checks ASCII digits and uint64 overflow before
each multiply/add. Leading zeros remain accepted. Name checks use unsigned bytes
with whitespace classification and reject empty/duplicate/NUL names. Empty contigs
retain zero line widths, require an in-file-or-EOF offset, and cannot enter positive
query division.

For a nonempty contig, `bases > 0`, `bytes >= bases`, and the first offset is inside
the opened FASTA. Let `q=floor((size-1)/bases)`, `r=(size-1) mod bases`, and
`available=file_bytes-1-offset`. Checking `q <= available/bytes` proves the product
`q*bytes` is representable and at most `available`; the subsequent subtraction and
tail check are safe. Because `bytes >= bases`, addressed positions are monotone,
so the final base bound also bounds every earlier base. Large unused strides for
a single addressed row remain valid and are explicitly tested.

`sequence()` rejects nonpositive positions before subtracting one, then rejects
`zero > size` or `n > size-zero` before reserving memory or performing line
arithmetic. A positive request to an empty contig is refused before division;
zero-length requests safely return empty. For accepted requests, every batched
address is bounded by the verified final address. No additional in-scope numeric,
range or valid-layout defect was found.

The only changed production file is `src/preservation_core.h`. All R APIs,
generated wrappers, projection code and reference writer formats are unchanged.
The complete core suffix from `inline Normalized normalize(` to end is byte-exact
to the delivered baseline, SHA-256
`59ab9d077a3d2fe3e7809bcdd566505ea7a30073120bd7c172c7338560020d67`.
Normalizer/parser/projection scientific versions remain unchanged; historical
measurements retain their historical source/DLL pins.

## Public behavior and meaningful RED

Independent R tests use format 2 and typed format 3, multiallelic SNPs, ID/contig
aliases, swapped measured alleles, an indel and a proven ALT–ALT contrast. They
compare complete old/new in-memory objects at batch sizes 1 and 4, all 19 normalized
columns, exact source tokens, literal effects, all nine ordinary edition files and
complete edition manifests after removing only creation time/output location.
Reference bundle bytes remain unchanged during these valid comparisons.

Additional deliberately re-pinned malformed FAIs refuse through ordinary opening,
explicit FASTA relocation and rebinding of a previously cached handle. A changed
on-disk FAI invalidates an existing handle even with `verify=FALSE`; explicit
rebinding then reads the updated trusted manifest and refuses the native unsafe
layout. Maximum unused strides and empty EOF offsets still produce exact baseline
normalized results.

The independent baseline program safely probes constructor acceptance and queries
only oracle-safe layouts. Against `392d826`, it records **10,503 constructor
mismatches** and **769 oversized-query exceptions**; the frozen candidate has zero
of either. Accepted malformed baseline layouts are never queried. This RED is
retained as `bounds-baseline-final.log`, separate from the implementation's own
original/native/public RED evidence.

An initial independent R expectation called a proven ALT–ALT result `ambiguous`;
the existing exact status is `alt_alt_contrast`. Complete old/new results already
matched. Four assertions against that overly broad literal were corrected without
changing source; `public-probes-discovery.log` and its summary retain the initial
result. The final independent public run is 86/86 with no errors or warnings.

On this macOS toolchain LeakSanitizer's `detect_leaks=1` is unsupported, and a
sandboxed `/usr/bin/time -l` could not read `kern.clockrate`. Those initial host-tool
diagnostics are preserved. Final ASan/UBSan programs were rerun directly with leak
detection disabled, both returning exit 0 with no ASan/UBSan finding. No claim of a
LeakSanitizer pass or timing/memory benchmark is made.

## Retained qualification and source integrity

All **868** candidate archive files and **132** exported baseline R/native files
were checked unchanged after testing. The candidate author checkout remained at
the frozen revision, clean, with `git diff --check` passing. No author source,
installed library, live full reference, preparation process, user repository or
cloud resource was modified.

The supplied author report hash was independently verified as
`645c12923a510d96e91d75c0b494f6b84746fa792580134c08d1d6fe9ce6e6ea`,
and author evidence index as
`df7aa173f5db492ae6a37026701604c533371fab3113c8f4ac127a67ba758571`.
All **756** indexed evidence files matched their sizes/SHA-256 values. The retained
qualification includes per-format 16,128 SNP and 79,904 haplotype rows, 297 explicit
scientific/resource cases, and byte-exact 96,032-row complete native output parity
to the previously reviewed scientific7 corpus. Those larger mathematical corpora
were checked as pinned retained evidence rather than rerun in this review.

The coordinator separately installed the exact frozen archive and ran the full
package suite. I independently checked its CSV/summary/guard and **868 matching
source pins**: **257 blocks / 2,570 passing assertions, zero failures/errors**, two
known warnings and two missing-reference skips. Its installed DLL is pinned as
`6f23605e06ac0b35fd40268b08190fd55e9a09e2fccded2f8cd7b629a3cd49eb`.
This is verified coordinator evidence, not an independent rerun of that larger
suite. Actual native Windows/macOS/Linux CI was pending at review completion.

## Reproduction and evidence

Evidence root: `results/mss-fai-safety-independent-review-stage8/`. Exact sources,
inputs, results and pins are in `bounds.cpp`, `make-bounds.py`,
`bounds-cases.tsv`, `bounds-input-summary.json`, `bounds-asan-final.log`,
`bounds-baseline-final.log`, `compatibility.cpp`, `compatibility-asan.log`,
`public-probes.R`, `run-public-probes.R`, `public-probes-final.log`,
`public-probe-summary.tsv`, portable/native logs, discovery/host-tool diagnostics,
`scientific-scope-verification.json`, `full-installed-verification.json`,
`final-source-integrity.json`, `final-summary.json` and `evidence-index.json`.

From that directory, representative independent commands are:

```sh
python3 make-bounds.py
clang++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Icandidate/src bounds.cpp -o bounds-asan
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 ./bounds-asan bounds-case-final bounds-cases.tsv
clang++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Icandidate/src compatibility.cpp -o compatibility-asan
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 ./compatibility-asan compatibility-final
Rscript run-public-probes.R .
```

The public harness loads local Rcpp/testthat from `/private/tmp/harvest-Rlib`;
adjust only that environment location when reproducing elsewhere. Portable suite
commands use `MSS_TEST_RLIB=/private/tmp/harvest-Rlib` with the archive's actual
SourceCpp/registered entrypoints. The native suite uses
`candidate/tests/native/run.sh`. Sanitizer options apply to this tested platform;
they are not new runtime dependencies of MSS.
