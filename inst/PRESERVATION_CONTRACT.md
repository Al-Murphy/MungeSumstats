# Preservation normalizer 7: scientific contract

This recipe describes `harmonise_sumstats()` and `harmonise_file()` using
`mss-preservation-7`, reference format `mss-preservation-2`, and the unchanged
`mss-projection-2` and parser identities. The legacy `format_sumstats()` API is
separate. Exact source, reference and algorithm versions belong in every receipt.

Every source row, identifier and original measurement token is retained. A
nullable normalized edition expresses what the supplied evidence supports;
uncertain normalized identity, signed effect, EAF and direction flags remain
missing. Retaining a source observation does not assert a verified canonical
association contrast.

## Measurement roles and source position

The generic `other` and `effect_allele` columns describe measured alleles. Neither
column declares FASTA REF or source strand. Only explicitly signed additive beta
or log-odds effects can be negated; raw odds ratios require caller conversion.
A frequency is transformed only when explicitly declared EAF. A swapped EAF is
complemented only when the catalog context establishes a biallelic event.

Generic complex DNA contrasts require a reported source position. This includes
indels and multi-base substitutions, even if one rsID identifies only one
canonical event. The same measured strings can admit different source anchors
and even opposite signs for one canonical event. With no position, distinct valid
complex DNA alleles return `source_coordinates_missing`, with all normalized
identity/effect/EAF fields withheld. Unsupported markers and identical alleles
retain their existing unsupported behavior; specialist entities remain preserved.
Single-base SNP ID inference remains supported: two distinct concrete single
bases cannot trim or shift, so canonical event positions exhaust their REF–ALT
anchors. Trusted source REF/ALT/strand declarations retain their stricter separate
coordinate requirements.

## Mapping pool and eligible ambiguity-proof universe

1. The original canonical mapping pool uses the existing coordinate and input-ID
   lookups and their intersection/conflict rules. Selection and
   `candidate_count` refer to this pool. The count is not a proof-graph size or
   a catalog-wide absence statement; a missing-position complex row can have
   positive candidate count while remaining unresolved.
2. Ambiguity evidence starts from the **full** coordinate and input-ID pools,
   before their intersection. Every current primary ID in their
   `Candidate.identifiers` sets participates; no first ID is selected.
3. Current-ID lookups expand to a fixed point. Exact canonical event coordinate
   lookups can find additional current IDs for the same event, including for
   ID-only queries. Only an already eligible exact canonical event bridges at
   this step. Unrelated event IDs returned from such a lookup do not become
   eligible merely because they overlap its coordinate. Initial full coordinate
   lookup candidates remain eligible as described in step 2.
4. Known source chromosome filters unrelated candidates before expansion. Input-ID
   and already loaded canonical coordinate lookups are reused. IDs, events and
   original context locations are deduplicated; cycles terminate.
5. Original S contexts are inspected separately by chromosome, position and REF
   span. Split rows at the same position and REF may combine ALT lists. Different
   REF spans or placements never pool their ALT strings. Original REF must be
   concrete DNA and validate against FASTA.

The completed eligible proof universe is a **rejection** mechanism. A distinct
REF–ALT interpretation `(event, swap, complement)` or a genuine paired ALT–ALT
contrast can make an otherwise selected mapping ambiguous. Proof-only mappings
never replace the main selected mapping or change its candidate-count meaning.
Known-position REF–ALT hypotheses use that exact source anchor; missing-position
single-base hypotheses use each closed event's canonical position.

This universe does not establish closure over every possible overlapping
haplotype in the catalog. Overlapping contexts with unlinked IDs and no eligible
coordinate/event bridge are outside it. It also cannot establish completeness of
a supplied reference or of missing study metadata. These limits must accompany
claims of harmonization. Original source provenance and explicit trustworthy
REF/ALT/strand contracts are preferred evidence.

## Equivalent ALT contrasts retain their backgrounds

An ALT–ALT witness requires both alleles in one FASTA-valid original REF context.
Joint common-prefix/suffix trimming retains the shared alternate background,
including flanks that differ from FASTA. Repeat-equivalent anchors compare both
original ALT events under one common inferred source REF span. No Cartesian
product of all ALT pairs is enumerated.

A shared alternate background is not erased merely because trimming yields the
same substitution and numeric beta sign as a reference-anchored event. For
example, original REF AC with ALT GC,GT compares GC against GT on a shared G
background. A shortened C/T observation can still represent that ALT–ALT contrast.
Without a resolving trusted source contract, the normalizer withholds its
reference-anchored effect rather than assigning the FASTA background silently.

## Proven ambiguity and incomplete context are different outcomes

`ambiguous` means that a competing REF–ALT interpretation or an exact paired
ALT–ALT witness rejects the selected generic mapping. `alt_alt_contrast` denotes
a proven ALT–ALT witness without a selected canonical REF–ALT mapping.
`source_coordinates_missing` defers a valid complex generic observation whose
source position is absent.

`source_allele_context_required` is a separate conservative outcome. It means
that complete background information is required to exclude a possible alternate
contrast; it does **not** assert that a necessary-only pair is the true contrast:

- For complex measured DNA alleles, a compatible pair has two distinct concrete
  ALTs whose absolute length difference equals that of the observed alleles.
  Joint padding/trimming and repeat reanchoring preserve this difference. If no
  such pair exists, this particular ALT–ALT explanation is impossible. A
  compatible pair is a necessary condition, not a complete haplotype proof.
- For single-base observations, a compatible concrete pair has equal-length ALT
  sequences differing at exactly one base with the observed base pair or its
  reverse complement. A shared-prefix/suffix join finds these without assigning
  a literal genome anchor or positive inferred REF span. Shared background
  insertions/deletions can change both quantities. A context that jointly reduces
  original REF and all concrete ALTs to a pure SNP has one proved SNP anchor;
  known source position different from that anchor excludes this context. The
  common flanks must match original REF, so this exception does not erase an
  alternate inserted background.
- An eligible FASTA-valid context with at least two distinct non-REF ALT entries
  and any nonconcrete ALT leaves genuinely unknown sequence or length. A missing
  `.` entry counts as unknown when the same context has another ALT. Identical
  ALT strings and strings equal to original REF do not supply a second allele.

Exact witnesses still take priority over conservative deferral. All three
unresolved outcomes withhold normalized identity, effect, EAF and direction
flags. They retain the source row and its measurements. Trusted complete source
REF/ALT/strand contracts keep their separate existing behavior. Ordinary
noncompeting three-allelic SNPs remain supported when their ALT base pair cannot
match the observation in either orientation.

For example, REF CAA / ALT GA,GAAA,CAAAA on a CAAAA… reference can make observed
A/AAA at a shifted repeat position look like a canonical insertion of two As.
The same strings can compare deletion of one A versus insertion of one A on a
shared G background. The fallback defers this generic interpretation rather than
claiming its FASTA-background effect. Explicit source REF A / ALT AAA metadata
would establish a different, trusted reference-anchored contrast and can resolve.

This is an explicit progressive-harmonization boundary. It retains variants and
studies with incomplete context and prevents them from silently supplying model
training labels. It does not make a completeness claim about specialist,
background/CIGAR or unlinked-overlap reconstruction.

## Resource limits and portability

`candidate_limit` continues to cap raw records per individual disk lookup. In
addition, the ambiguity-proof closure admits at most that many unique current
IDs, canonical events, original context locations, and examined S REF contexts.
Additional closure work is bounded by **16 × candidate_limit** units, counting
raw records before lookup deduplication, ALT candidates, candidate signature
positions and normalization attempts. The bound is derived explicitly from the
configured limit; it is not an unbounded alias graph or ALT-pair search.

Exhaustion raises `Reference context closure limit exceeded`; it never silently
skips remaining evidence to publish a selected effect. Existing raw lookup
exhaustion keeps its candidate-limit error. No partially completed file edition
is published. References and original acquisition objects remain available for
retry with a reviewed higher resource limit.

Closure stores current-ID/event/context locators and consumes owned lookup
batches, without retaining copies of complete ALT-payload vectors. Signature
views live only while their owning allele strings remain alive. Work remains
sensitive to allele-string length, original input-line size, and reference repeat
length; candidate counts are not a byte-level or wall-time guarantee. The normal
bounded parser/batch/reference caches remain in effect. Operator-scale runs
should additionally enforce process-tree memory, disk and time bounds.

The required implementation remains portable CPU C++17/R. It does not require
GPU processing, WSL, containers, an external bcftools installation or a new
reference format.

### Opt-in typed reference data format 3

`prepare_reference(..., reference_format="mss-preservation-3")` is an opt-in
storage prototype. The default remains `mss-preservation-2`. The normalizer is
still `mss-preservation-7`, the delimited parser is `mss-delimited-1`, and the
projection algorithm is `mss-projection-2`; data-format selection does not add
scientific evidence or relax any nullable result contract.

The same immutable directory files are used. `core-version` is exactly
`mss-preservation-3` followed by LF. Each LF-terminated `records` row is either:

```
A <chr> <canonical_pos> <REF> <ALT> <current_id> <original_pos>
S <chr> <original_pos> <original_REF> <complete_ALT_list>
```

Fields are separated by one TAB. A has seven fields and S has five, including
the tag. Every parsed source-row occurrence has its own S offset, independently
of whether any sequence ALT is supported. Every successfully normalized ALT
occurrence has its own A offset. Identical physical row contents are not globally
deduplicated. Sorting removes only duplicate `(lookup_key, same_offset)` entries,
so candidate-limit checks still count raw lookup occurrences before semantic
merging. Primary and alias I keys are unchanged; A stores the current primary ID,
and aliases remain opaque UTF-8 keys. S never supplies identifiers to the
scientific closure. The complete original six-column TSV remains byte-exact in
`source.tsv`; record REF/ALT text retains the existing uppercase representation.
Rejected/nonsequence/empty ALT contexts remain in S and source, and rejected ALT
occurrences remain in `rejected.tsv`. Signed nonpositive original positions may
occur in S, as in format 2; A positions must both be positive.

C/I reads decode only A. The original source position locates S when paired site
evidence is required. S reconstructs its own original chromosome, position and
REF with an empty candidate ALT before semantic merging, preserving distinct REF
spans. No first-A record doubles as S, and no C/I candidate eagerly owns or unions
the complete site ALT list. Readers support formats 2 and 3; earlier format-2-only
readers reject 3. Format 3 validates complete selected typed records and key rows:
field/tag/key type and coordinate agreement, concrete distinct canonical A DNA,
canonical integer spelling and range, UTF-8 without NUL/TAB/LF framing bytes,
in-range offsets at physical row boundaries, and complete terminating LF,
including lines across cache pages. It does not certify unread unrelated records
or the global sort order through a query; immutable artifact checksums remain the
whole-file integrity check. Lookup alias membership cannot be reconstructed from
A alone; the checksummed key index and original source are authoritative.

Format-3 manifests and output receipts carry `record_layout="typed-A7-S5"`
(reference manifests) and `reference_context_order="semantic-s-context-v1"`.
After raw lookup caps and same-span merging, S contexts are examined in semantic
`(chromosome, original_position, original_REF, merged_ALT_list)` order. The first
three fields uniquely identify a merged context; thus no physical offset, source
row ID, omitted S identifier, or ALT-list tie selects context order. ALT token
occurrence/order within an unmerged source context remains unchanged, including
repeated raw tokens charged to the scientific work budget. This versioned resource
order may produce a nullable uncertainty result where physical-offset format 2
refused an exhausted proof budget, or a refusal where format 2 encountered an early
witness. Such bounded resource outcomes are classified explicitly; exhausted proof
never publishes an arbitrary canonical effect. Default format-2 reader behavior
and writer payload bytes remain unchanged. Adequate-budget scientific rows must
agree across both formats.

Projection supports all four source/target format combinations. Existing 2-to-2
receipts retain their fields; any combination containing 3 adds explicit
`source_reference_format` and `target_reference_format`. `reference_format` is
the common format when equal and `mixed` otherwise. Neither coordinate build nor
normalization/projection algorithm is inferred from the data format. Format 3
requires independent scientific/storage review, native Windows CI and complete
reference capacity admission before default or full production adoption.
