#!/usr/bin/env bash
set -euo pipefail
[[ $(uname -s) == Linux ]]
packet=$(cd "$(dirname "$0")" && pwd -P)
[[ "$packet" == /workspace/mss-batch-rank-cost-diagnostic-20261122 ]]
prior=/workspace/mss-batch-rank-benchmark-source-20261122
old=/workspace/mss-key-page-8k-experiment-20261122
inputs=/workspace/mss-task4/inputs
results=/workspace/mss-batch-rank-cost-results-20261122
export LC_ALL=C
cd "$packet"
mkdir "$results"
sha256sum -c SOURCE.sha256 > "$results/source-before.txt"
(cd "$prior" && sha256sum -c SOURCE.sha256) >> "$results/source-before.txt"
c++ --version > "$results/compiler.txt"
for role in baseline candidate; do
    extra=(); if [[ $role == candidate ]]; then extra=(-DMSS_BATCH_RANK_BENCH); fi
    c++ -std=c++17 -O2 -g -DMSS_LOOKUP_COUNTERS "${extra[@]}" \
        "-I$prior/$role" "-I$old" driver.cpp -o "$results/$role"
done
(cd "$results" && sha256sum baseline candidate) > "$results/binaries.sha256"
for order in matched shuffled; do
    for role in baseline candidate; do
        prefix="$results/$order-$role"
        python3 "$old/meter.py" "$prefix.process.json" "$results/$role" panel \
            "$inputs/reference/matched-index" "$inputs/reference/GRCh37.ensembl104.chr22.fa" "$order" \
            < "$inputs/mss/matched-1M.tsv" 2> "$prefix.timing.json" |
            python3 "$old/digest.py" 1000000 > "$prefix.digest.json"
        cmp "$prefix.digest.json" "$old/expected/panel.json"
    done
done
sha256sum -c SOURCE.sha256 > "$results/source-after.txt"
(cd "$prior" && sha256sum -c SOURCE.sha256) >> "$results/source-after.txt"
