# Paired original-site ALT diagnostics — MSS implementation 5 candidate

Source change: `ba91c5598083c3511a91da81f3ce3b7bc66f7d54`, branch
`fix/paired-site-contexts`, isolated clone
`/private/tmp/gwas-harvest-stage2/mss-site-context-v5`, parent frozen MSS
`5cf1e23c907ad9c923141faf0728964951be7c66`.

This is a **scientific uncertainty-diagnostic fix**, separate from any reference
packing change. The reference format/default writer remains
`mss-preservation-2`; no compact writer, schema 3, compression or source deletion
was introduced. The frozen source, installed libraries and active full-reference
preparations were not modified. The candidate is ready for independent review;
full installed-package validation and native Windows CI remain later gates.

## Concrete issue and final behavior

Before the change, a failed allele match could pool ALT strings from different
original positions/REF spans and compare that union against the first merged
candidate's REF. Textual record-offset ordering could therefore change
`alt_alt_contrast` into `allele_mismatch` without changing the logical reference
sites. Two ordinary frozen-format-2 preparations from the same rows in different
orders reproduce the issue. Both outcomes were unresolved in the demonstrated
case; no resolved coordinates, beta direction or EAF correction was gained.

The failed-match diagnostic now considers each original S context separately.
An ALT pair must occur together in that original position/REF span, the REF must
be concrete DNA and validate against FASTA, and any known input chromosome and
position must agree with that context. Split rows at the **same** position and
REF span remain merged by the existing S lookup, retaining real multiallelic
site context. Unsupported-only, symbolic and historical-alias context remains
in the reference; no source row or reference record is dropped.

Only the `hits.empty()` uncertainty path changed in the C++ normalizer. Matching,
allele flipping, EAF complementation, canonical normalization, candidate lookup
and reference writing retain their existing algorithms. The normalizer
implementation identity is now `mss-preservation-5` in in-memory provenance,
file metadata/manifest and projection normalization receipts. The projection
algorithm identity remains `mss-projection-2`; the index format remains 2.

Reviewed candidate-source SHA256s:

| File | SHA256 |
| --- | --- |
| `src/preservation_core.h` | `1409b4fa2b1b3a156f9abd1d7b8ba00c4076eccecabae2fa576e22e44659fd85` |
| `R/preservation.R` | `14bf6e6df65640609baed03ea87c066abd1b84fd90439f73e3eb4eaf7d46d14d` |
| `R/projection.R` | `ef0bd49102df563f6f9c2f8284a233b5a8af3014a6a493e48ba2185e57c5acad` |

## Observed RED and targeted GREEN

The new public-R tests were run against the frozen installed FASTA-offset build
before implementation. **Seven functional assertions failed for the intended
reason**: false pooled `alt_alt_contrast` labels. These include three of six
source-order permutations, different REF spans at one coordinate, borrowing a
pair from a different original position/chromosome, and unsupported REF context.
The genuine paired-site positives and resolved-direction controls already passed.

After implementation, the complete portable standalone preservation runner
passed **37 blocks / 463 assertions**, including eight new gold blocks / 38
assertions. That runner also exercises existing reanchored indel context,
source REF/ALT contracts, aliases, symbolic context, quoted streaming,
projection, contig aliases and FASTA layouts. It compiles/loads the candidate
native source and R wrappers without loading the whole Bioconductor package.

The new gold cases preserve genuine same-context ALT-ALT pairs, including
symbolic third-allele context; retain source input exactly; leave unresolved
effects/frequencies nullable; and retain match/swap/complement/swapped-complement
beta and biallelic EAF behavior. `git diff --check` passed.

Raw receipts: `results/paired-site-context-v5/red.log` and
`standalone-green.log`. This is not a claim that the full installed package suite
or Windows CI has run for implementation 5.

## Independent old/new resolved-field comparison

A bounded native functional matrix uses the same frozen-format-2 index for both
compiled readers, with legacy small-reference cases plus reanchored indels,
multiple IDs/aliases, multiple contexts, unsupported/symbolic alleles, unknown
coordinates, both generic and declared REF/ALT contracts, three frequency kinds,
and deliberately small candidate caps. It is not a speed or full-genome benchmark.

| Measurement | Result |
| --- | ---: |
| Total functional queries | 190,392 |
| Generic queries | 175,392 |
| Declared-contract queries | 12,000 |
| Deliberately bounded-cap queries | 3,000 |
| Resolved rows with **all public fields byte-exact** | 4,548 |
| Complete result rows byte-exact | 190,061 |
| Identical candidate-cap errors | 815 |
| `alt_alt_contrast` → `allele_mismatch`, only status changed | 327 |
| Old cap error → unresolved `allele_mismatch` | 4 |
| New errors / new or lost resolved mappings | 0 / 0 |

For the 327 successful diagnostic changes, every non-status field is exact,
including candidate count and direction flags. The four avoided errors all use
`candidate_limit=1`, unknown position and a known chromosome different from the
ID candidate's chromosome: the new diagnostic skips the unrelated chromosome's
S lookup instead of hitting its two-context cap. They remain unresolved, with
no beta/EAF or coordinate assignment. No filtering/mapping/sign gain is claimed.

The initial comparison deliberately treated any error/success difference as a
failure and exposed those four cases. They were inspected and classified; the
final classifier rejects new errors, any resolved-field change, any successful
non-status change, and any status transition outside the diagnostic pair above.

An actual tiny writer parity check also built the same reference with old/new
native default writers: all six native payloads (`core-version`, FAI, keys,
offsets, records and rejected rows) were byte-identical. R's source-TSV copy and
reference-manifest preparation code was not changed.

Evidence under `results/paired-site-context-v5/`:

- `compare_driver.cpp`, `compare_resolved.py`, `classify_existing.py`;
- `resolved-comparison/result.json`, with input/output SHA256s;
- `resolved-comparison/format2-writer-parity.json`;
- exact tiny reference, generated queries and both full outputs remain locally
  in that results directory; they are not copied into the source repository.

## Remaining gates

Obtain independent scientific/code review before adoption. Then build/install a
separate candidate library, run the full installed-package suite and actual
Windows/macOS/Linux native CI, and update Harvest's explicitly pinned normalizer
identity/config only when that code is reviewed. Existing implementation-4
benchmark/CI evidence keeps its original identity.

The optional compact-reference discussion remains in
`reports/mss-reference-capacity-feasibility-stage5.md`. Even if this correction
removes incidental-order dependence in the audited diagnostic path, compact
reader reconstruction and all outputs still need their own proof. No compact
format is approved or delivered by this change.
