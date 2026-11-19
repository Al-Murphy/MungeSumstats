#!/usr/bin/env bash
# $1 is an empty, dedicated scratch root containing only the staged inputs,
# source archives, and sanitized admission note. Run as a child of guard.py.
set -eu -o pipefail
ops=$(cd -- "$(dirname -- "$0")" && pwd)
bench=$(dirname "$ops")
repo=$(cd "$bench/../../.." && pwd)
scratch=${1:?absolute scratch path}
deps=${2:?immutable image shared R library}
candidate=${3:?candidate commit}
image=${4:?image digest}
base=a86e8c94d753630eff8c6a695ed2f240d1123d0d
archives="$scratch/archives"
result="$scratch/results"
mkdir "$result"
run_logged() {
    local log=$1
    shift
    "$@" >"$result/$log.log" 2>&1
}
sha_for() {
    awk -v name="$1.tar" '$2==name {found++; digest=$1} END {if(found!=1) exit 1; print digest}' "$archives/SHA256SUMS"
}
base_sha=$(sha_for "$base")
candidate_sha=$(sha_for "$candidate")
run_logged build python3 -B "$ops/build.py" "$result/build" --candidate "$candidate" \
    --archives "$archives" --baseline-sha256 "$base_sha" --candidate-sha256 "$candidate_sha" \
    --dependencies "$deps" --image-digest "$image"
installations="$result/build/installations.json"
run_logged smoke-prep python3 -B "$ops/assets.py" smoke "$installations" "$result/small" \
    --admission-note "$scratch/admission.txt"
run_logged smoke python3 -B "$bench/run.py" run "$result/small/smoke.json" \
    "$result/smoke" --suite diagnostic
run_logged functional bash "$ops/functional.sh" "$result/build" "$deps" "$result/small" "$result/functional"
run_logged full-prep python3 -B "$ops/assets.py" full "$installations" "$result/assets" \
    --inputs "$scratch/inputs" --admission-note "$scratch/admission.txt"
matrix_rc=0
run_logged matrix python3 -B "$bench/run.py" run "$result/assets/matrix.json" "$result/matrix" || matrix_rc=$?
echo "$matrix_rc" >"$result/matrix-exit.txt"
case "$matrix_rc" in 0|2|3) ;; *) exit "$matrix_rc";; esac
for name in distinct-shuffled coordinate-only id-only negative dense mixed hard-cases declared; do
    for role in baseline candidate; do
        f2="$result/matrix/$name-f2/pair-01-$role/edition/normalized.tsv"
        f3="$result/matrix/$name-f3/pair-01-$role/edition/normalized.tsv"
        run_logged "$name-$role-cross-format" cmp "$f2" "$f3"
    done
done
for scale in 1M 10M 40M; do
    run_logged "snp-oracle-$scale" python3 -B "$repo/inst/verification/fai-safety-stage8/whole-file/harness/verify.py" \
        --profile "$result/assets/snp-oracle.json" --stage "$scale" \
        --edition "mss7=$result/matrix/matched-$scale/pair-01-candidate/edition" \
        --receipt "$result/snp-oracle-$scale.json"
done
run_logged legacy python3 -B "$bench/run.py" run "$result/assets/legacy.json" "$result/legacy" --suite legacy
run_logged alloc python3 -B "$bench/run.py" run "$result/assets/diagnostic.json" "$result/alloc" \
    --suite diagnostic --profile alloc --allocation-threshold 65536
run_logged cpu python3 -B "$bench/run.py" run "$result/assets/diagnostic.json" "$result/cpu" \
    --suite diagnostic --profile cpu
run_logged lookup-pressure "$result/functional/candidate-lookup" pressure "$result/functional/counter-fixture"
exit "$matrix_rc"
