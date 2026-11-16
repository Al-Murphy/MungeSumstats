# Frozen MSS7 whole-file benchmark

Preparation is separate from execution. `selfcheck.py` exercises the independent streaming oracle and strict workload boundary. `bench.py preflight --name <name>` seals this harness, verifies every reference payload against the RDS receipt, checks source/library pins and prior input MD5s, and writes an immutable profile. Root supplies a separately pinned SCIENCE_GREEN approval for that profile/source/DLL/full installed suite.

Run only an authorized stage with `bench.py run-stage --profile <profile> --approval <approval> --stage 1M`. Larger stages require the exact matched 1M pair to have passed plus explicit approval. Each new profile/stage is append-only; old failed evidence is never overwritten.

The baseline is frozen implementation4, not legacy `format_sumstats`; previous legacy measurements stay historical. The candidate is implementation7 with reference-format2 and projection not requested. Only this scientific version's correctness overhead is compared here. Report package setup, reference verification, processing and full process separately. One R/data.table thread; a fresh process/handle per case; uncontrolled OS cache on a shared host.

Inputs repeat 437,458 real chr22 unique biallelic non-palindromic SNPs with synthetic exact beta/EAF. The 40M workload is 40M rows, not unique variants. An independent FASTA-validated source oracle examines all 19 normalized fields of every row and both exact source copies. Broader reference contexts are deliberately refused by this workload oracle, not approximated; indels/multiallelic sign and identity correctness are covered by the separate reviewed scientific matrices.

Every owned subprocess tree is monitored for RSS, macOS physical footprint or Linux swap, with an 8GiB limit, 160GiB free-space floor, finite time and owned-output-growth limits. No unrelated process is signaled. External `time` peaks are also checked. Setup/output/verifier resource use and source retention remain separate denominators. Monitoring cannot guarantee that an unobserved sub-poll transient never occurred.

Version2 closes the independent review of v1: exact anchored scientific GREEN source verdict, installed-suite source revision binding, undeclared/missing row-width refusal and strict unique locus/primary-ID denominator. Version1, its seal/profile and adverse evidence remain unchanged; no v1 whole-file benchmark was executed.
