The matched benchmark uses the legacy parent b4a2d21 and the installed current
package. Run `benchmark-matched.R mode input reference output metrics cache threads`
under external `/usr/bin/time -l` (macOS) or `-v` (Linux). `mode` is `new` or
`legacy`; `cache` is `cold` or `warm`. Each process records timing phases and input
identity. Reference-cold means no live reference handle at entry, **not** a flushed
filesystem cache. Warm means a verified reusable native handle or loaded legacy
reference data.table, primed by one complete identical file pass. External process
time/RSS includes package startup, reference initialization and priming; use the
recorded processing phase for the core-plus-file costs.

For an exact scientific comparison, supply biallelic, non-palindromic SNP records
with known coordinates, rsIDs, alleles, signed BETA and EAF. The reference TSV and
FASTA are exactly those used to prepare the new disk reference. The legacy
reference acquisition function alone is replaced with a data.table containing
those exact reference records; installed parent QC functions are unchanged.
This is an adapted legacy reference baseline, **not** an unmodified end-to-end
SNPlocs retrieval benchmark. The adapter's loading and representation construction
are timed separately. The unmodified installed package test suite is also run.

`on_ref_genome=FALSE` avoids the legacy unknown-ID imputation branch for this
all-known panel. `infer_eff_direction=FALSE` fixes the declared A2 effect-allele
meaning; `check_dups=FALSE` permits repetitions of the real panel for scale;
`frq_is_maf=FALSE` keeps declared EAF semantics. Coordinate sorting is disabled in
both paths. Other QC, allele matching/flipping and biallelic filtering execute.
Matched input repeats are disclosed rather than represented as independent
biological variants. Added functionality is measured on a separate challenge
panel, never credited as matched SNP speedup.

`run-matched-benchmarks.py` supervises the full sequential matrix with immutable
runner snapshots, per-case OS time/RSS and disk receipts, a32GiB process-tree
memory guard, and complete output correctness scans. A guard termination is
reported as censored evidence and is never called a natural OOM. Warm process
time/RSS includes priming; phase metrics separate priming from the timed pass.
The verifier checks every retained row, normalized REF/ALT, signed BETA and EAF.

On macOS the guard must use `max(RSS, phys_footprint)`, including compressed and
private accounting; RSS alone is inadequate. On Linux it uses RSS plus VmSwap.
The monitor polls at100ms and fails closed if a live process cannot be measured.
An actual16MiB/112MiB child allocation test proves under/over64MiB behavior;
the complete supervisor also has a64MiB startup-censorship receipt. These are
sampled safety thresholds, not guarantees against allocation-granularity overshoot.
