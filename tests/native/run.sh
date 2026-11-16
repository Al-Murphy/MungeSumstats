#!/bin/sh
set -eu
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
"${CXX:-c++}" -std=c++17 -Wall -Wextra -pedantic -O2 "$(dirname "$0")/test_core.cpp" -o "$work/core-test"
"$work/core-test" "$work/fixtures"
"${CXX:-c++}" -std=c++17 -Wall -Wextra -pedantic -O2 "$(dirname "$0")/test_blocks.cpp" -o "$work/blocks-test"
"$work/blocks-test" "$work"

"${CXX:-c++}" -std=c++17 -Wall -Wextra -pedantic -O2 "$(dirname "$0")/projection.cpp" -o "$work/projection-test"
"$work/projection-test" "$work/projection.chain"

"${CXX:-c++}" -std=c++17 -Wall -Wextra -pedantic -O2 "$(dirname "$0")/test_format3.cpp" -o "$work/format3-test"
"$work/format3-test" "$work"
