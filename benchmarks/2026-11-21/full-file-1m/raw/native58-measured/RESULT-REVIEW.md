# Native58 measured result review

**The measured invocation exited zero in 6.04 s. All 58 methods, 14 phase outcomes, 21 case records and scoped process checks passed review.** This run supplies the independently observed foreground duration missing from fix2. The root observation is a draft awaiting operator review and recording. The matrix consumer binding and interpreter alias handling remain pending.

The root operator ran `/usr/bin/time -l /bin/zsh /private/tmp/mss-fullfile-20261121/native58-measured/run-native.sh` and observed foreground exit 0. The initial tool result was `17bdc2`, session `52892`, and the terminal result was `a3ef39`. I verified the retained artifacts and timing receipt. This review did not launch, rerun, or signal any tested process, and wrote only this review and the draft observation. The original Stage185 and fix2 evidence roots remain separate from this run. [Measured source review](REVIEW.md) · [Prior fix2 review](../native58-fix2/RESULT-REVIEW.md)

## Timing, resources and roster

- The operator stdout is empty and the status has exact bytes `0\n`. The stderr contains one timing row with **6.04 s real, 0.54 s user and 0.28 s system**. The timed utility is the two-argument zsh launcher call. The measured duration is positive and below 180 s. The suite self elapsed is 5.82144487509504 s, within that measured duration and its own 130 s budget. The suite self peak RSS is 43,941,888 B, below 1.5 GiB. The observation uses that suite RSS field. [Operator timing](../native58-operator/foreground.stderr) · [Suite summary](green/suite.stdout)
- The pinned source declarations, loader and legacy method roster agree on 58 unique methods across nine classes. The verbose log contains exactly those 58 ordered `ok` lines and the expected final summary. The suite JSON reports exit 0 with zero failures, errors or skips. The suite wrapper stdout is empty and the saved suite status is `0\n`. [Verbose results](green/suite.stderr)
- The run recorded 297,106,084 KiB free, above the 178,257,920 KiB floor. The inspector's 147,770 logical B and 348,160 allocated B for `green/cases` match an independent metadata sum and both fall below 32 MiB. The largest retained regular file under `green` is 8,843 B, below 32,768 B. [Capacity receipt](green/admission.txt) · [Closure report](green/owned-process-closure.json)
- All 12 manifest entries verify. The entry, before and after receipts match all 12 successful pin checks. The four-file source-freeze and source-recheck match the manifest and current bytes. The wrapper hash is `8b710fd6795e1016784fa72a7e25b9f23fdb132afddd002022bbb6cfef42d6a6`, the manifest hash is `90dbd6c02ea3280f678686daf22c6258ab13ef13623787d95c6dcaa0d0f01a87`, and the timing utility hash is `435da47a8eb633f605aa895ed41e270c3907085b97b5b5c67aa9df1011c41895`. [Manifest](native-pins.sha256) · [Source freeze](green/source-freeze.sha256)

## Workflow and scoped evidence

All four workflow steps match their recorded command arrays and distinct original child PIDs. Each reports ordinary `waitid` code 1, status 0, raw `wait4` status 0 and decoded exit 0. The recorded ownership command uses the Homebrew `opt` alias. It and the manifest's `/opt/homebrew/bin/python3` both resolve to the physical `/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15` and match SHA-256 `6cc907e4a022fefc843bcd93954eaf5a2cbfc77a095640bfadeeed96154190c5`. [Workflow receipt](green/workflow-status.json)

| Step | PID | Recorded command |
| --- | ---: | --- |
| `pins_before` | 48886 | `/usr/bin/shasum -a 256 -c native-pins.sha256` |
| `suite` | 48887 | `/bin/zsh run-tests.sh green` |
| `ownership` | 49098 | `/opt/homebrew/opt/python@3.15/bin/python3.15 -IBu native-verify.py --inspect` |
| `pins_after` | 49100 | `/usr/bin/shasum -a 256 -c native-pins.sha256` |

I matched each inspector phase to its planned, running and terminal journal entries, including the exact leader PID, group and session, closure/reap, sampled RSS and expected diagnostics.

| Phase | Leader | Raw exit | Primary | Closure |
| --- | ---: | ---: | --- | --- |
| `fast` | 48971 | 0 | None | confirmed |
| `nonzero` | 48974 | 7 | `phase_exit` | confirmed |
| `allowed-two` | 48977 | 2 | None | confirmed |
| `signal` | 48980 | -15 | `phase_exit` | confirmed |
| `foreground` | 48983 | 0 | None | confirmed |
| `leftover` | 48999 | 0 | `leftover_descendants` | confirmed |
| `rss` | 49025 | -15 | `rss` | confirmed |
| `SIGINT` | 49030 | -15 | `shutdown` | confirmed |
| `SIGTERM` | 49034 | -15 | `shutdown` | confirmed |
| `SIGHUP` | 49038 | -15 | `shutdown` | confirmed |
| `term-aware` | 49042 | 17 | `wall` | confirmed |
| `lost` | 49088 | -15 | `ps_timeout` | unresolved |
| `transition` | 49093 | 23 | None | confirmed |
| `descendant-transition` | 48903 | 31 | None | confirmed |

- The term-aware child used 1.0 s wall and sample intervals, took one live sample, and recorded primary `wall` with an empty cleanup list. It received exactly one sent group SIGTERM after the deadline. It was observed as a matching zombie and reaped at exit 17. The signal record contains no SIGKILL. Stdout is exactly `ready\nterm observed\n`. [Term-aware journal](green/cases/test_term_aware_child_closes_before_escalation/term-aware.jsonl) · [Term-aware stdout](green/cases/test_term_aware_child_closes_before_escalation/term-aware.stdout)
- The intentional `lost` phase records unresolved phase closure and the exact `ps_timeout` primary and sole cleanup entry with detail `synthetic query timeout`. Its recorded PID and group are covered by final scoped absence. No required phase records `ps_helper_cleanup`. The separate patched observer regression tests intentional synthetic helper-cleanup reporting. [Lost journal](green/cases/test_uncertain_closure_uses_only_direct_child_and_blocks_pending/lost.jsonl) · [Patched observer regression](test_guard.py#L676)
- The decoy records agree on PID/PGID/SID 48998, alive immediately after guard cleanup. Leftover PID 49001 belonged to group 48999 and appears in the group's verified SIGKILL target list after SIGTERM. The descendant readiness and wait files confirm the parent reaped PID 48905 at exit 19. The supplement records one injected transient row, a fresh native observation and the outer leader reap `[48903, 31]`. Both unfinished descendant temp paths are absent. The direct metrics children report the required exits and complete Darwin wait4 RSS, PID 49089 at exit 0 with 30,261,248 B, and PID 49090 at exit -15 with 13,910,016 B. [Leftover journal](green/cases/test_leftover_descendant_escalates_kill_and_decoy_is_left_alive/leftover.jsonl) · [Descendant wait](green/cases/test_parent_reaps_descendant_before_fresh_native_query/parent-wait.json) · [Descendant supplement](green/cases/test_parent_reaps_descendant_before_fresh_native_query/injection.json)

I independently derived 20 PIDs and 14 groups from the phase leaders, sampled, terminal and signal members, and the decoy, leftover descendant, two direct metrics children and reaped descendant. The resulting exact command matched both the inspector and its kernel wait record.

```text
/bin/ps -x -p 48903,48905,48971,48974,48977,48980,48983,48985,48998,48999,49001,49025,49030,49034,49038,49042,49088,49089,49090,49093 -g 48903,48971,48974,48977,48980,48983,48998,48999,49025,49030,49034,49038,49042,49088,49093 -o pid=,pgid=,ruid=,uid=,state=,rss=
```

Scoped query PID 49099 had ordinary `waitid` code 1 and status 1, raw `wait4` status 256 and decoded exit 1. Both streams were empty, `known_absent` and `complete` were true, and `incomplete` was empty. [Scoped closure receipt](green/owned-process-closure.json)

## Draft status and remaining matrix work

The [draft observation](DRAFTroot-observation.json) has the exact Stage185 v1 schema, the measured two-argument launcher command, operator-observed exit zero, 6.04 s elapsed, the suite peak RSS and available KiB, and the actual 11 source, 12 runtime artifact and 21 case path/size/SHA-256 pins from this run. Draft SHA-256 is `0e5b7e6dd33608980315cdd872239e39d446fc6b6a9c4423afd0d2b68b5a9cfb`. The draft awaits a separate operator review and independent recording outside the packet and matrix run trees. Its current in-packet path is for review.

The legacy matrix consumer is hard-coded to the original packet root, command and source digests. Its ownership check accepts the literal launcher's Python path or the full Cellar path, while the authentic wait evidence records the verified Homebrew `opt` alias. Acceptance of this draft's measured run therefore still requires reviewed binding and alias handling. The pending plan has null functional, native and primary proof slots, and both matrix run roots were absent during this review. No matrix activation or scientific run is part of this result. [Stage185 root observation schema](/Users/bschilder/code/gwas-harvest/ops/mss-observer/v2/README.md#L145) · [Legacy proof consumer](/private/tmp/mss-performance-20261119/tools/darwin-matrix-stage185/native-proof.py#L60) · [Pending matrix plan](/private/tmp/mss-performance-20261119/tools/darwin-matrix-stage185/plan.pending.json)

## Operator and supplemental pins

The exact v1 object has an exhaustive set of fields and 12 native artifact keys. The three separate operator receipts are pinned here under `/private/tmp/mss-fullfile-20261121/native58-operator`.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `foreground.stdout` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `foreground.stderr` | 777 | `c25f5d09b49755639cd5a78ae7fc1313325d736158577179510c7403d8d5b36d` |
| `foreground.status` | 2 | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

Two supplemental focused records under `green/cases` were also verified.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `term-aware.stdout` | 20 | `eda7ac330720dc2795d89cbbb0907df676cc5c5f0dfb89a36a45656f32b8a379` |
| `injection.json` | 213 | `3d9fd9ebb604619ff2814b2fda0c782bbfb1eba5ec8279263f33b2908ed3f9db` |
