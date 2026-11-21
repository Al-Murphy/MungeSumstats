# Prepare the matched 1M diagnostic

[prepare-1m.py](prepare-1m.py) has SHA-256 `1b33e6ce562c6b9a294987f458e38abac066124a48e031cb57a5c1c182cd9b63` (7,716 B). It reads the actual current Darwin installations and the accepted rebuilt format-2 reference, verifies the relevant frozen harness and oracle sources, and produces three new metadata files under `fullcall-1m-01`.

| Bound input | SHA-256 |
| --- | --- |
| `installations.json` | `8b5b4d9a104f34da2053e21c99f6e235021088ac695741704f815284ba27578b` |
| `dependency-inventory.json` | `293025cbce6ea392e6b6ecbda1787b90b904abd5a2f30958a9301505fcf46a31` |
| `reference-accepted.json` | `dce55f634f56040ab5072955351ce05b0191ffdcad16b9cda0a9487ebdbc36aa` |
| 712-byte rebuilt manifest | `d15ad358b438fef485ddf0980043556a43f3e4d6b72faca35f3429ce4919ab93` |
| Recovered FASTA | `f2e991e0816c2d7a19f9252ba314ea2717ed1ca73563dcd6417d514fb6c85e53` |
| Ordered reference panel | `d6fae3acdc98e4d06f2e5eb01492360aeb874fe5b9f8763914cc39d9c40db9d1` |
| `matched-1M.tsv` | `d110828a8b17bf88309db042e7073faa267a0b63a3550d14a87f7761074c91f2` |

The adapter rechecks the common and framework package identities, both installed role payloads and build receipts, all seven stable reference members and the pinned new manifest. The operator's manifest acceptance binds both successful public open receipts and their stdout. The frozen ordered-panel verifier compares the entire input to the expected alternating records and confirms 1,000,000 rows from 437,458 panel events. The frozen runner checks schema, asset hashes, installed provenance and cardinality using a deep copy.

After scoped smoke qualification and with a current operator note describing observed load, serial measurement, and the uncontrolled desktop file cache, the publication command is

```sh
TOP=/private/tmp/mss-fullfile-20261121
/opt/homebrew/opt/python@3.15/bin/python3.15 -I -B -S \
  "$TOP/runtime/prepare-1m.py" --admission-note "$TOP/fullcall-1m-01/admission-1m.txt"
```

The new `diagnostic-1m.json` specifies one `matched-1M` case, three alternating pairs, explicit verified FASTA binding, EAF, beta, tab, and the frozen batch settings. It retains the actual note, 186 GiB continuing floor and 2 GiB estimated output. It uses `target: false`. Run it with the frozen runner's unprofiled diagnostic mode. The independent `snp-profile.json` binds the ordered reference panel, FASTA, FAI, build and input for the later external oracle. `matched-1M.keys.json` supplies the input-bound cardinality receipt. All three paths use exclusive writes.

A source-only `--check` run exited 0 and printed 1,000,000 rows, 437,458 events and three diagnostic pairs. The temporary cardinality scratch was removed, and an independent check found all three publication paths absent. The publication path, R benchmarking, and independent SNP output verification remain unexecuted in this task.
