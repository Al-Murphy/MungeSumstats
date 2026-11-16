# Final source MSS scale results

**Timing comparison correction:** read the [authoritative erratum](stage17-final-mss-scale-timing-erratum.md) before the frozen report. Final versus historical legacy 10M external wall is 229.370 versus 142.250 seconds; processing is 224.500 versus 137.107 seconds. The reports formerly mixed these scopes. Receipt values and the frozen 47-file archive/tar are unchanged.

Final production `03d11e92b417add0183dea4fb77716f5502d5735`, measured DLL `6f23605e06ac0b35fd40268b08190fd55e9a09e2fccded2f8cd7b629a3cd49eb`: 1M/10M/40M rows all pass every 19-field check and both source-copy checks. Processing is 24.444/224.500/811.171s; peak R RSS is 1.518/1.858/2.254GiB. Repeated 437,458-record chr22 SNP workload, single runs/shared host/uncontrolled page cache; no legacy runtime speedup, unique40M or whole-genome equivalence claim.

See [complete table, exact evidence and limits](stage17-final-mss-scale/REPORT.md), [release pins](stage17-final-mss-scale-release.json) and [compact archive](stage17-final-mss-scale.tar.gz). Existing counter40 and prior94 evidence are unchanged. No new benchmark, full count/profile/admission, index, provider action or production supervision qualification is authorized. Full reference preparation/count/index and full corpus remain outstanding. Root owns the combined main-plan endappend.
