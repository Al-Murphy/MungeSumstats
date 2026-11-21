# Eliminating full-width copies used to count allele flips

The median paired speedup was **1.42×** across four allele-flip calls on separate 100,000-row synthetic GWAS inputs with a supplied reference. The revised function counts matching row positions, reuses the counts and keeps the original selectors for allele swaps, effect updates and imputation. The four calls took **129–148 ms** in the control and **91–111 ms** in the treatment. Separate profiling recorded **30.42% fewer allocated bytes**.

| Pair | Process order | Control ms | Treatment ms | Control / treatment |
| --- | --- | ---: | ---: | ---: |
| 1 | c1, t1 | 129 | 91 | 1.418 |
| 2 | t2, c2 | 148 | 111 | 1.333 |
| 3 | c3, t3 | 148 | 94 | 1.574 |

Each worker processed the four cells in the order below. Timing covered one function call per cell. Fixture preparation, garbage collection and result checks took place outside the timer. Narrow inputs had 11 columns and wide inputs had 27. Sparse and dense inputs flipped 1,000 and 75,000 rows respectively. Every cell retained 1,000 mismatches, enabled the imputation flag and used all six numeric effect columns with mismatch logging and dropping disabled.

| Width / density | Control ms, pairs 1/2/3 | Treatment ms, pairs 1/2/3 | Control allocated B | Treatment allocated B | Reduction |
| --- | --- | --- | ---: | ---: | ---: |
| Narrow / sparse | 16 / 19 / 19 | 13 / 16 / 15 | 22,323,200 | 20,438,224 | 8.444% |
| Wide / dense | 68 / 78 / 78 | 45 / 55 / 42 | 198,370,912 | 129,073,072 | 34.933% |
| Narrow / dense | 32 / 35 / 36 | 22 / 27 / 26 | 102,288,944 | 71,395,968 | 30.202% |
| Wide / sparse | 13 / 16 / 15 | 11 / 13 / 11 | 20,474,912 | 18,073,072 | 11.731% |
| Total | 129 / 148 / 148 | 91 / 111 / 94 | 343,457,968 | 238,980,336 | 30.419% |

The allocated-byte totals matched across all three repetitions of each role. BSD time measured whole-process peak RSS of 341.8–448.7 MB for timing and 363.3–456.8 MB for profiling. The treatment's timing peak exceeded the control's in pairs 1 and 2. These source-level results apply to the tested Apple Silicon R process, supplied reference and synthetic inputs. Installed-package execution, production data, full-pipeline throughput and Windows performance remain unmeasured.

## Correctness and runtime

The original source came from commit `737ad24784ee9a2445a3eaf4f8262c87561a72f7`. The control passed 25 exact oracle and literal cases. Its work assertion then counted four full-width row subsets used to obtain counts, giving the expected exit status 1. The treatment passed the same 25 cases and produced zero such subsets. The trace observed one full-width missing-j control, zero explicit-j controls and zero narrow controls. It restored the data.table namespace, S3 method and options and passed the traced parity check. Unit comparisons cover return values and attributes, caller mutations, condition classes/messages/calls and log bytes. Two partial-flip fixtures exercise a coercion warning from an integer `tmp` column and an assignment error from an existing factor `flipped` flag. All four 100,000-row fixtures also passed exact control/treatment comparison before measurement.

The historical runtime used arm64 macOS, R 4.5.3, data.table 1.18.6.1, one data.table thread, one BLAS/OpenMP thread and disabled JIT. The unit, verification and 12 measurement processes each exited zero with matching entry/exit checks of source, framework, tools and data.table files. Full SHA-256 manifests are in `pins/`.

| File | SHA-256 |
| --- | --- |
| Original R/check_allele_flip.R | `50161aa8c60ba2c48ddb0d6dc95114487da6325f9218ba0decda0f5dd0c007b2` |
| Revised R/check_allele_flip.R | `0859ac4809ccac9bb2b14f2971571e02e02dfda0cd5d511da27de92dc649f6e9` |
| test_flip_indices.R | `b485f5522efb7c9e7171c0358a67958f4c51b9b542a5ea67b6abd225176d60aa` |
| benchmark_flip_indices.R | `fce09a27d0c8e86f8d94100b949b49d9df1315ad57df5e446e3ae8f28d1067fc` |
| data.table shared library | `0f578ffd675064a727f8b0a2b9cf3b3891f0e493354676a6fe209f936ac017c6` |
| Historical launcher | `83965c1ba2ddf883b124820eeaa6355b0cc2aec6756c14d3c8b0907b1b2ba189` |

`baseline/`, `unit/` and `verify/` contain the worker receipts. The baseline and treatment used the same 25-case test harness. The `time/` and `alloc/` directories contain raw stdout, stderr and exit status per trial. Time-worker stdout reports R user, system and elapsed seconds. Allocation-worker stdout reports counted Rprofmem events and allocated bytes. The archive contains all 24 original profiles. For example, the following reads the wide/dense control profile from this directory and reproduces its 845 counted events and 198,370,912 bytes.

```sh
tar -xOzf allocation-profiles.tar.gz c1/out/wide.dense.Rprofmem |
    awk '/^[0-9]+ / {events++; bytes += $1} END {printf "%d %.0f\n", events, bytes}'
```

## Reproducing the workers

Run from the repository root at the treatment revision on a machine with the matching R framework, data.table installation and macOS time/timeout utilities. Set `MSS_R_FRAMEWORK` and `MSS_DATA_TABLE` to the corresponding physical directories. The commands use a fresh working directory, reconstruct the four original R sources from Git and relocate the recorded launcher and manifests. The launcher verifies the source hashes and runtime files before and after each invocation.

```sh
set -eu
export MSS_R_FRAMEWORK="${MSS_R_FRAMEWORK:-/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources}"
export MSS_DATA_TABLE="${MSS_DATA_TABLE:-$HOME/Library/R/arm64/4.5/library/data.table}"
evidence="$PWD/benchmarks/2026-11-21/flip-indices"
work=$(mktemp -d /tmp/mss-flip.XXXXXXXX)
mkdir -p "$work/oracle/R" "$work/r-component-red-v1" "$work/r-component-red-v3"
ln -s "$PWD" "$work/r-worktree"
for file in check_allele_flip.R messager.R get_unique_name_log_file.R write_sumstats.R; do
    git show "737ad24784ee9a2445a3eaf4f8262c87561a72f7:R/$file" > "$work/oracle/R/$file"
done
cp "$evidence/pins/"*.before.sha256 "$work/r-component-red-v1/"
cp "$evidence/pins/tools-pins.sha256" "$work/r-component-red-v3/"
cp "$evidence/run-r-component.sh" "$work/"
relocate() {
    sed -e "s|/private/tmp/mss-focus-20261120|$work|g" \
        -e "s|/private/tmp/mss-performance-20261119/parser-precedence-stage176|$work/oracle|g"
}
relocate < "$evidence/pins/r-green-source-pins.sha256" > "$work/r-green-source-pins.sha256"
relocate < "$work/run-r-component.sh" |
    sed -e 's|^RFW=.*|RFW="$MSS_R_FRAMEWORK"|' -e 's|^DT=.*|DT="$MSS_DATA_TABLE"|' > "$work/replay.sh"
baseline_status=0
sh "$work/replay.sh" "$work/baseline" unit control || baseline_status=$?
test "$baseline_status" -eq 1
grep -F "GWAS row subsets with logging and drop disabled: 4 (expected zero)" "$work/baseline/worker.stderr"
sh "$work/replay.sh" "$work/unit" unit treatment
sh "$work/replay.sh" "$work/verify" verify control
for phase in time alloc; do
    mkdir "$work/$phase"
    for trial in c1 t1 t2 c2 c3 t3; do
        case "$trial" in c*) role=control ;; t*) role=treatment ;; esac
        sh "$work/replay.sh" "$work/$phase/$trial" "$phase" "$role"
    done
done
```
