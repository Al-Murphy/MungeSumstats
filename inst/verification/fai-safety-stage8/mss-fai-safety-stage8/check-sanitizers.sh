#!/bin/sh
set -eu
source_root=/private/tmp/gwas-harvest-stage2/mss-fai-safety-stage8
work=/private/tmp/gwas-harvest-stage2/results/mss-fai-safety-stage8/sanitizers
mkdir "$work"
clang++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer "$source_root/tests/native/test_fai.cpp" -o "$work/test-fai"
"$work/test-fai" "$work"
