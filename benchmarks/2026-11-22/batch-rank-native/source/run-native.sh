#!/usr/bin/env bash
set -euo pipefail
if [[ $(uname -s) != Linux ]]; then echo "Requires Linux" >&2; exit 2; fi
packet=$(cd "$(dirname "$0")" && pwd -P)
[[ "$packet" == /workspace/mss-batch-rank-benchmark-source-20261122 ]]
old=/workspace/mss-key-page-8k-experiment-20261122
inputs=/workspace/mss-task4/inputs
results=/workspace/mss-batch-rank-benchmark-results-20261122
export LC_ALL=C
cd "$packet"
mkdir "$results"
sha256sum -c SOURCE.sha256 > "$results/source-before.txt"
(cd "$old" && sha256sum -c SOURCE.sha256) >> "$results/source-before.txt"
(cd "$inputs" && sha256sum -c "$old/panel-inputs.sha256") >> "$results/source-before.txt"
c++ --version > "$results/compiler.txt"
flags=(-std=c++17 -O2 -g)
for role in baseline candidate; do
    extra=(); if [[ $role == candidate ]]; then extra=(-DMSS_BATCH_RANK_BENCH); fi
    c++ "${flags[@]}" "${extra[@]}" "-I$role" "-I$old" driver.cpp -o "$results/$role"
    c++ "${flags[@]}" -DMSS_LOOKUP_COUNTERS "-I$role" "$old/lookup_parity.cpp" -o "$results/parity-$role"
done
(cd "$results" && sha256sum baseline candidate parity-baseline parity-candidate > binaries.sha256)
fixture="$results/mixed"
"$results/baseline" prepare-fixture "$fixture"
for role in baseline candidate; do
    "$results/parity-$role" parity "$results/parity-fixture-$role" > "$results/lookup-parity-$role.txt"
    cmp "$old"/expected/lookup-parity.txt "$results/lookup-parity-$role.txt"
    for order in matched shuffled; do
        for format in 2 3; do
            ref="$fixture/format$format"
            prefix="$results/fixture-format$format-$order-$role"
            "$results/$role" fixture "$ref" "$fixture/tiny.fa" "$order" 2>"$prefix.timing.json" |
                python3 "$old/digest.py" 12 > "$prefix.json"
            cmp "$prefix.json" "$old/expected/fixture-format$format.json"
            for mode in short empty; do
                prefix="$results/$mode-format$format-$order-$role"
                "$results/$role" "$mode" "$ref" "$fixture/tiny.fa" "$order" < "$old/controls/$mode.tsv" \
                    > "$prefix.rows.tsv" 2> "$prefix.timing.json"
                if [[ $mode == short ]]; then
                    cmp "$prefix.rows.tsv" "$old/expected/short.rows.tsv"
                else
                    test ! -s "$prefix.rows.tsv"
                fi
            done
        done
    done
done
python3 summarize.py "$results" controls > "$results/controls.json"
for order in matched shuffled; do
    for pair in 1 2 3; do
        roles=(baseline candidate)
        if [[ $pair == 2 ]]; then roles=(candidate baseline); fi
        for role in "${roles[@]}"; do
            prefix="$results/panel-$order-$pair-$role"
            python3 "$old/meter.py" "$prefix.process.json" "$results/$role" panel \
                "$inputs/reference/matched-index" "$inputs/reference/GRCh37.ensembl104.chr22.fa" "$order" \
                < "$inputs/mss/matched-1M.tsv" 2> "$prefix.timing.json" |
                python3 "$old/digest.py" 1000000 > "$prefix.json"
            cmp "$prefix.json" "$old/expected/panel.json"
        done
    done
done
python3 summarize.py "$results" paired > "$results/paired.json"

sha256sum -c SOURCE.sha256 > "$results/source-after.txt"
(cd "$old" && sha256sum -c SOURCE.sha256) >> "$results/source-after.txt"
(cd "$inputs" && sha256sum -c "$old/panel-inputs.sha256") >> "$results/source-after.txt"
