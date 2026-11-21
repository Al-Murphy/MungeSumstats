# Installed-suite stop handling re-review

**Specification compliance — pass. Code quality — pass.** The original Important finding is addressed. Open findings in this fix are 0 Critical, 0 Important and 0 Minor.

The operator records the phase outcome and verification state, then checks the frozen guard before its first post-run evidence read. An expired deadline or pending shutdown at phase return goes directly to the funded failure-result path. During verification, the operator checks the guard as it traverses and reads installed files, reads staged tests, resolves reference links, verifies runtime controls and parses receipt rows. The accepted writer applies its ordinary stop checks when publishing the post-run inventory. A later stop can retain an already completed inventory and record incomplete verification. A closed test failure with time remaining still retains the available evidence and exits nonzero. [Phase return and post-run verification](key-geometry-installed-suite.py#L326) · [Installed inventory](key-geometry-installed-suite.py#L115) · [Receipt parsing](key-geometry-installed-suite.py#L186)

The revised receipt parser consumes the condition rows once and retains the previous warning and error calculations. The reviewed controls exercise both stop reasons at phase return, shutdown during verification and publication, and an ordinary closed failure. The author reported a zero exit and recovered the prior 269 cases and 3,077 passing assertions. The controls simulate the child outcome and use the actual frozen stop checks and metadata writer. This review did not execute R or repeat payload hashing. The fresh suite result remains pending. [Controls](key-geometry-installed-suite-fix1-controls.py) · [Reported command and results](KEY-GEOMETRY-INSTALLED-SUITE-FIX1-REPORT.md)

| Verified artifact | SHA-256 |
| --- | --- |
| Fixed operator | `d06a1286dfcd183ea00aa1a5f2d8c73061b44bf9cbdf9a35912444aab8a75a91` |
| Fix diff | `72e65c80e67421aff93b741d1da083fa42693c4edee89fca58856acc7aab82b5` |
| Controls | `00f1d63e10f5f5cf84db9c19ca59598bf9aad93b20e41b25ee9bc82cba719601` |
| Fix report | `01e81f432582b48538e023db0d14e688111fd343acd7d38a6a102be0bbcc42df` |
| Original source review | `87ea651e681a177ee7ba21f55f3eda449c522e1c552812b8068d8e188bd1a82f` |
