#!/bin/sh
# Local component launcher; fixed runtime, sources and retained receipts.
set -eu
umask 077
ulimit -c 0
ulimit -f 4096
RROOT=$1
MODE=$2
ROLE=$3
case "$MODE" in unit|verify|time|alloc) ;; *) exit 2 ;; esac
case "$ROLE" in control|treatment) ;; *) exit 2 ;; esac
WT=/private/tmp/mss-focus-20261120/r-worktree
CONTROL=/private/tmp/mss-performance-20261119/parser-precedence-stage176/R
RFW=/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources
DT=/Users/bschilder/Library/R/arm64/4.5/library/data.table
PINS=/private/tmp/mss-focus-20261120/r-green-source-pins.sha256
BASE=/private/tmp/mss-focus-20261120/r-component-red-v1
TOOLS=/private/tmp/mss-focus-20261120/r-component-red-v3/tools-pins.sha256
mkdir -m 700 "$RROOT"
mkdir -m 700 "$RROOT/tmp" "$RROOT/out" "$RROOT/lib"
ln -s "$DT" "$RROOT/lib/data.table"
cd "$RROOT"
cp "$PINS" source-freeze.sha256
(cd "$RFW"; /usr/bin/shasum -a 256 -c "$BASE/framework.before.sha256") > framework-entry-pins.txt
(cd "$DT"; /usr/bin/shasum -a 256 -c "$BASE/datatable.before.sha256") > datatable-entry-pins.txt
(cd "$RFW/library/tools"; /usr/bin/shasum -a 256 -c "$TOOLS") > tools-entry-pins.txt
/usr/bin/shasum -a 256 -c source-freeze.sha256 > source-entry-pins.txt
case "$ROLE" in control) SOURCE="$CONTROL/check_allele_flip.R" ;; treatment) SOURCE="$WT/R/check_allele_flip.R" ;; esac
if [ "$MODE" = unit ]; then
    set -- "$WT/tests/standalone/test_flip_indices.R" "$CONTROL/check_allele_flip.R" "$SOURCE" "$CONTROL" "$RROOT/lib" 1.18.6.1 "$DT/libs/data_table.so" "$RFW" 4.5.3 "$DT" "$RROOT/out"
else
    OTHER=unused
    if [ "$MODE" = verify ]; then OTHER="$WT/R/check_allele_flip.R"; fi
    set -- "$WT/tests/standalone/benchmark_flip_indices.R" "$MODE" "$WT/tests/standalone/test_flip_indices.R" "$SOURCE" "$OTHER" "$RROOT/out" "$RROOT/lib" 1.18.6.1 "$DT/libs/data_table.so" "$RFW" 4.5.3 "$DT"
fi
if /usr/bin/time -l /usr/bin/timeout -s TERM -k 1s 45s /usr/bin/env -i PATH=/usr/bin:/bin LC_ALL=C TZ=UTC TMPDIR="$RROOT/tmp" R_LIBS_USER=NULL R_LIBS_SITE=NULL R_ENABLE_JIT=0 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "$RFW/bin/Rscript" --vanilla --default-packages=methods,utils "$@" > worker.stdout 2> worker.stderr
then
    worker_status=0
else
    worker_status=$?
fi
echo "$worker_status" > worker.status
(cd "$RFW"; /usr/bin/shasum -a 256 -c "$BASE/framework.before.sha256") > framework-exit-pins.txt
(cd "$DT"; /usr/bin/shasum -a 256 -c "$BASE/datatable.before.sha256") > datatable-exit-pins.txt
(cd "$RFW/library/tools"; /usr/bin/shasum -a 256 -c "$TOOLS") > tools-exit-pins.txt
/usr/bin/shasum -a 256 -c source-freeze.sha256 > source-exit-pins.txt
/usr/bin/awk '/maximum resident set size/ {if ($1 >= 1073741824) exit 1}' worker.stderr
/usr/bin/du -sk . > retained-usage.txt
exit "$worker_status"
