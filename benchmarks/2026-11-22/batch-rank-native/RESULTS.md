# Native batch-rank experiment: retain the existing lookup

The prototype passed the scientific checks but was slower in all six paired comparisons. It is not selected for Rcpp/R integration.

| Median paired candidate change | Matched input order | Shuffled input order |
| --- | ---: | ---: |
| Native preparation plus harmonization | +7.4144% | +5.5069% |
| Whole-worker wall time | +5.1357% | +4.3098% |
| Whole-worker CPU | +6.4482% | +4.7735% |
| Peak worker RSS | −0.0121% | −2.0245% |

Positive time differences mean slower. Three fresh serial AB/BA/AB pairs ran for each order, twelve workers total, using the same 1M repeated-SNP input and reference. Both roles use the selected16KiB×512 key cache. The candidate clock includes key gathering, estimation, prefetch sorting/search, original-order harmonization, directory clear and destruction of caller key storage. No R API timing or broader genome-wide throughput is inferred.

Selected MSS revision is25c26fe2d86441360541f4b505bf546405bb14f2; baseline header5881233fbad917d4475f4801497c90a720428c41ba5afb8a7d38c234da60f155 and experimental header9f71fa92eaa8e6b8d249ba58d7863630a7f6958e491ea413a48d521151932225. Neither header changed during the run. All twelve ordered13-field1M digests match the frozen independent gold, including effect/frequency binary64 bits. Both reference formats passed the mixed, short and empty controls, and both native parity binaries matched historical output. All ten candidate prefetch attempts admitted in each panel process. The maximum core arena receipt is5,945,737B; the9,600,000B caller estimate is separate and is not a strict heap bound.

The original operator stopped before compilation because its1MiB source-stage assumption was smaller than MFS directory charges. A fresh reviewed operator reused the exact staged payload and kept stage allocation within8MiB and the unchanged256MiB phase limit. Its native shell finished successfully, but the guard subsequently exited1 because a tuple in six expected execution-order fields was compared directly to lists decoded from JSON. All92 result metadata files returned unchanged. Root independently reproduced this exact mismatch and verified exact equality after JSON projection, then recomputed the paired arithmetic independently with Decimal.

A separate verification-only phase completed at07:19:21UTC, exit0, without repeating any scientific worker or compiler. It verified all92 metadata hashes, actual native exit0 and unsignaled quiescent drain from the original guard journal, complete gold/order controls, post-source/input/binary hashes and original storage accounting. Proof SHA256:dc7beb39b434ded686030652cc98ad994293c751fa43993fe3aa32c969bb884e. The original failed terminal remains unchanged. Four unpersisted guard-local fields cannot be recovered: original compiler-before SHA, original input-before identities, exact monitor-return RSS peak and final storage sample preceding the tuple failure. The separate proof labels these limitations and records current checks; it does not manufacture the original success proof.

The measured outcome rejects this prototype on this workload. The cause of its extra CPU cost requires separate instrumented attribution; the timing alone does not establish that sorting is responsible. The selected R-only cleanup and all broader full-matrix, real-source and Windows gates remain separate.
