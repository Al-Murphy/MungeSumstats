#!/bin/sh
set -eu
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
"${CXX:-c++}" -std=c++17 -Wall -Wextra -pedantic -O2 "$(dirname "$0")/test_core.cpp" -o "$work/core-test"
"$work/core-test" "$work/fixtures"
