# Native count-only reference refinement

Implements the separately requested exact refinement from
`reports/full-reference-format3-capacity-stage7.md`. This is an operator diagnostic,
not an MSS API/default change or permission to build a production index.

1. Differential tests build only small synthetic format2/3 fixtures using the
   frozen production writer, then independently measure records/rejected/keys and
   offsets. Include duplicate keys, aliases, multiallelic/unsupported contexts,
   indels, invalid positions, Unicode, CRLF and repeated historical rows.
2. A no-output C++ pass reuses the unchanged pinned normalizer. Count exact record
   offsets, raw sort bytes and within-record unique keys; distinct record offsets
   prove cross-record deduplication is impossible. Checked arithmetic refuses
   overflow. Report maxima and all raw/rejected occurrences without filtering.
3. A trusted coordinator verifies exact source/header/FASTA hashes before and
   after, owns resource guards and computes retained/merge envelopes only with an
   explicit bound for small files/manifest. Small-file allowance, allocation,
   simultaneous work and free-space reserve remain separate.
4. Independent review precedes use on any sealed full input. Neither active
   preparation nor a full native index is read/modified/started by development.

Frozen production header SHA256:
`4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc`.
No production header, normalizer, R API, wrappers or installed DLL is edited.
