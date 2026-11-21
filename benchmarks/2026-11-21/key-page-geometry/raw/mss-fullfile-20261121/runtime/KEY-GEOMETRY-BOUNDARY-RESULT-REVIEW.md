# Typed key boundary independent result review

**Spec compliance passes. Evidence quality passes. Critical 0, Important 0, Minor 0.** Both guarded builds passed the three format-3 boundary cases, legacy parity, and legacy pressure checks. This closes the native boundary gate before the public-API trial. [Boundary brief](KEY-GEOMETRY-BOUNDARY-BRIEF.md) · [Accepted source review](KEY-GEOMETRY-BOUNDARY-REVIEW.md)

## Run evidence

- Root identified launch chunk `fcfcbe` and observed session `11263` exit `0` in terminal chunk `46abae`. I verified the root observation SHA-256 `c3252da11e0c4c6b53823284c821faa2a084f32a644fd7783564203dbc2e9507`, the result SHA-256 `e4fffc3e14de38c90a51ad0b2bed823868d51ad74ad6b53b34b4062d90f8a696`, and all six guard receipt hashes. Both compile, both parity and both pressure receipts have raw exit `0`, confirmed closure, reaped leaders, empty stderr and no signals, triggers or cleanup errors. The root observation records all six groups and sessions absent. [Root observation](/private/tmp/mss-key-geometry-experiment-20261121/root-observation-02-boundaries.json) · [Run result](/private/tmp/mss-key-geometry-experiment-20261121/run-02-boundaries/result.json)
- The compile receipts use the accepted sibling driver SHA-256 `b11f9d1897aacc32fc3899f546f27c6958a908a7df648f6bc1a5c24d37cf7888`, the same arm64 C++17 `-O2` flags and macOS 27.2 SDK, and `MSS_LOOKUP_COUNTERS` with assertions enabled. I verified the source review and both role-specific header pins. The parity and pressure plans invoke their corresponding binaries and separate fixture directories. Baseline and trial executable SHA-256 values are `018d65cce4caac3e42f4e5984332517ce0af18fd3755d1729783138fb1a294ba` and `e7d5daec74e7c0c30b218a346866f8a78f43b28a23df509c52ee08a2ff98f52f`. [Admission](/private/tmp/mss-key-geometry-experiment-20261121/run-02-boundaries/admission.json)

## Observed parity

I compared the two full 17,354-byte transcripts byte-for-byte. Both have SHA-256 `4401910e98e5bdffde63e4485ebaf83928cb0d678ab07492dbc8c47af3be09ee`. The first 16,585 bytes match the accepted run's format-2/3 parity transcript exactly. I checked every line of the appended transcript against the expected outcomes and inspected the three small generated key and offset layouts for each role. [Baseline transcript](/private/tmp/mss-key-geometry-experiment-20261121/run-02-boundaries/parity-baseline.stdout) · [Trial transcript](/private/tmp/mss-key-geometry-experiment-20261121/run-02-boundaries/parity-trial.stdout)

| Boundary case | Verified generated layout and observed result |
|---|---|
| Newline at page end | LF at byte 16,383 and second key row at 16,384. Both aliases returned one complete expected candidate and consumed one unit of context work. |
| Malformed typed offset | `I:r` points to the row carrying `01` at 16,384. The exact typed unsigned integer error was reported at both initial work values, with work still `0` and `1` respectively. |
| Split UTF-8 | `E9 9B AA` spans offsets 16,383–16,385 and the short key row starts at 16,389. Two alternating long-short pairs each returned the full expected candidate and consumed one unit per lookup. |

The retained format-3 fixtures use the expected typed `A` record and exact little-endian offsets in both roles. The pressure transcripts match each role's accepted prior transcript. Both formats report zero allocation calls for 1,024 warm-negative queries and complete the 3,000-distinct-query checks. [Baseline pressure](/private/tmp/mss-key-geometry-experiment-20261121/run-02-boundaries/pressure-baseline.stdout) · [Trial pressure](/private/tmp/mss-key-geometry-experiment-20261121/run-02-boundaries/pressure-trial.stdout)

This review inspected bounded retained artifacts without recompiling, rerunning the driver or hashing real reference payloads. The frozen counter review still has its original SHA-256 `d478ab3a248afb1f94cca78c913037735aadea29d8e6bd9c4fbbd539c4f425d0`. The boundary run qualifies correctness on these fixtures. Public-API timing, peak RSS and full 19-field parity remain separate gates. [Counter result review](KEY-GEOMETRY-RESULT-REVIEW.md)
