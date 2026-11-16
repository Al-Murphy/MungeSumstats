# Native FAI arithmetic safety

Close the inherited numeric-validation follow-up from the relocation review.
FAI field meanings follow https://www.htslib.org/doc/faidx.html . This is numeric
and bounds validation, not a claim to independently reconstruct/authenticate an
operator-rewritten FASTA/index bundle. All successful scientific outputs must stay
unchanged; invalid indexes fail before native arithmetic.

1. Add constructor-only negative native cases and valid LF/CRLF/no-final-newline,
   Unicode, empty-contig and leading-zero decimal controls; preserve baseline RED.
2. Parse decimal fields completely with checked uint64 overflow; reject empty or
   duplicate names, malformed rows and impossible widths/file bounds. Preserve
   the existing writer's zero-length contig representation without division.
3. Check requested ranges without unsigned addition overflow. Do not restrict
   terminators to LF/CRLF: the standard permits other consistent terminators.
4. Add R tests with deliberately re-pinned malformed FAI in both formats, plus
   exact valid output/receipt parity. Register new R tests in both portable paths.
5. Independent review, native/installed/oracle checks and source-pinned evidence
   before adoption. Existing throughput measurements keep their original pins.
