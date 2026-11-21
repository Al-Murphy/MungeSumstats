# Completed-scales verifier review

Reviewed source SHA256 `f6d90ac9d20af9a469f0f9098bc39ce220f8393bdb1a733cfdadd897d3edb14b`. The addition in `completed-scales-source.diff` reproduces the source byte for byte and its SHA256 is `4aa14bdab3c84432df047fd969487da3f68f9ea5d7874ef85ee23da4b1f7d2c0`. I read the brief, source, relevant frozen helper implementations, original small controls and author report. This was a source review. I did not execute the verifier, R or author controls, hash scientific payloads or traverse installed trees.

## Important finding

**Twelve retained worker logs are missing from the promised control inventory.** At `verify-completed-scales.py:164–165`, the explicit list contains only stdout and stderr. The frozen helper's `ROLE_LOGS` includes `r-heap.rds` at `warm-current.py:54`. The frozen `worker.R:241` writes this per-worker heap diagnostic even for profile `none`. Direct stat calls on the twelve expected original paths confirmed all twelve files exist and are 305–311 bytes. Their hashes are consequently absent from the closing checks at lines 279–280 and from `RESULT.json.original_small_controls` at line 308. A changed or missing retained heap log would go undetected. This affects evidence completeness, independently of the wait4 RSS measurement.

Replace the literal pair of filenames on line 164 with `warm.ROLE_LOGS`. The existing `controls.pin`, closing recheck and result serialization already supply the required behavior. The fix needs no new output, reads only the twelve small log files at each check, and leaves the frozen summary and thresholds intact.

## Spec and quality verdict

**Spec requires the control-inventory fix before invocation. Quality is otherwise acceptable for the scoped guarded experiment. No critical findings.**

- Lines 69–106 authenticate the disk-floor observation, planned/running records, supervisor exit 1, raw exit −15, confirmed closure and missing terminal record. The original post-pins and summary are expected to be absent. The status of the original matrix is preserved.
- Lines 206–245 bind the config, SNP profile, opening, source pins, role builds, live inventories and selected cardinalities. Each of the twelve observations comes from its original physical run directory, in the required alternating order. Frozen settings, namespace, runtime and process checks bind the role-specific installation, positive wait4 RSS and zero exit. Source and small controls are rechecked after the second live inventory.
- Lines 253–273 use genuine frozen byte comparison, digest and scan logic on all six editions at each scale. Frozen manifest verification runs the R manifest-only command, compares stable fields and requires a passing marker and successful process receipt. These new evidence files are written in fresh case directories under the operator's private output root. Source paths, original process outputs and original editions are read inputs.
- Lines 174–186 and 283–293 retain paired ratios, ratio of processing medians, whole-worker elapsed/CPU medians, RSS medians/maxima and relative reductions. Speed plus median RSS and the independent candidate 8 GiB absolute threshold have separate status fields. Verifier exit 0 reports successful execution of its own checks.
- Lines 189–203 and 294–311 bind each planned oracle window to the physical original pair-01 candidate edition, source/profile and full input/artifact digests. Both 19-field oracles and the root observer remain pending. The result limits claims to the selected scales and records warm ineligibility.
- The author report scopes the proposed invocation to its own 4 GiB RSS, 16 MiB metadata, 600 s wall and 186 GiB floor. There is no new 50 GB warm admission, full-matrix restart or runtime-success claim in the report.
