#!/usr/bin/env bash
# Run inside the owned Linux guard. No background processes.
set -eu -o pipefail
OPS=$(cd -- "$(dirname -- "$0")" && pwd)
build=${1:?build directory}
deps=${2:?shared dependency library}
science=${3:?small prepared inputs}
out=${4:?new functional output directory}
mkdir "$out"
export MAKEFLAGS=-j1 R_LIBS_SITE="$deps" R_MAKEVARS_USER=/dev/null
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
for role in baseline candidate; do
    source_dir="$build/$role/source"
    lib="$build/$role/library"
    export R_LIBS_USER="$lib:$deps"
    sh "$source_dir/tests/native/run.sh" >"$out/$role-native.log" 2>&1
    for suite in preservation registered; do
        Rscript --vanilla "$source_dir/tests/standalone/$suite.R" "$source_dir" \
            >"$out/$role-$suite.log" 2>&1
    done
    Rscript --vanilla -e \
      'a<-commandArgs(TRUE);.libPaths(c(a[2],a[3]),include.site=FALSE);testthat::test_dir(file.path(a[1],"tests/testthat"),package="MungeSumstats",load_package="installed",reporter="summary",stop_on_failure=TRUE)' \
      "$source_dir" "$lib" "$deps" >"$out/$role-testthat.log" 2>&1
    Rscript --vanilla "$OPS/failure_parity.R" "$lib" "$deps" "$science" "$out/$role-errors" \
        >"$out/$role-errors.log" 2>&1
    # One fixture source across both builds. Counters are candidate-only.
    flags=()
    if [[ "$role" == candidate ]]; then flags+=(-DMSS_LOOKUP_COUNTERS); fi
    c++ -std=c++17 -O1 "${flags[@]}" -I "$source_dir/src" \
        "$build/candidate/source/tests/standalone/lookup_parity.cpp" -o "$out/$role-lookup" \
        >"$out/$role-lookup-build.log" 2>&1
    "$out/$role-lookup" parity "$out/$role-lookup-fixture" \
        >"$out/$role-parity.stdout" 2>"$out/$role-parity.stderr"
done
cmp "$out/baseline-errors/errors.tsv" "$out/candidate-errors/errors.tsv"
cmp "$out/baseline-parity.stdout" "$out/candidate-parity.stdout"
