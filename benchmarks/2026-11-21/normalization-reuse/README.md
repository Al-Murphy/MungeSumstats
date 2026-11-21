# Row-local normalization reuse on a mixed native fixture

The treatment took **35.85–37.16 ms** for the 24 measured case loops, compared with **48.54–50.62 ms** for the control. The median paired speedup was **1.36×**. Each worker ran the same 12 synthetic rows through reference formats 2 and 3, with 256 consecutive calls per row and format. That gives 6,144 calls per process. These measurements describe warm, repeated calls to the C++ header on Darwin. Installed R, production-file throughput, Windows and full-reference workloads require separate measurement.

| Pair | Execution order | Control seconds | Treatment seconds | Control / treatment |
| --- | --- | ---: | ---: | ---: |
| 1 | c1, t1 | 0.048535584 | 0.036894122 | 1.315537 |
| 2 | t2, c2 | 0.050618332 | 0.035854459 | 1.411772 |
| 3 | c3, t3 | 0.050406958 | 0.037157291 | 1.356583 |

Each TSV gives the case, format, iteration count, wall seconds and last status. The table sums the 24 case times per process and excludes the two reference-open measurements. For the shifted repeat deletion, the ratio of median control time to median treatment time was 2.683× for format 2 and 2.690× for format 3. Generic indel and missing-coordinate cases showed small, inconsistent changes at submillisecond durations. The corresponding BSD time -l files include preparation, open, output and teardown. All six processes exited zero. Maximum resident set size ranged from 3,194,880 to 3,309,568 bytes, which does not establish a memory saving.

## Scientific and work checks

The control is commit `737ad24784ee9a2445a3eaf4f8262c87561a72f7`. The single treatment change retains complete source-anchor normalization values in four row-local slots with 256 bytes of owned character storage. The generic coordinate seed, hit matching and competing REF–ALT proof were counted separately from biallelic EAF normalization. The selected swapped biallelic SNP made 4/4/4/1 calls in control and 4/0/0/1 in treatment in both formats. The diagnostic compared full native result fields, lookup keys, returned candidate records, context debit order and exceptions before evaluating the work reduction.

| Check | Observed result | SHA-256 of output |
| --- | --- | --- |
| Both old/old roles and treatment, full diagnostic transcript | Byte-identical, 201,381 bytes each | d8cd96341d29e332d143a62d59950ae47fd5ee2c8c474e369dc2a9f6aa10a551 |
| Old/old work records | Each format saved 0 of 8 required, comparator exit 1 | c1c39639b30a8d614f22a214f48484ad1d966c64a0608fe899c3f55ad829384a |
| Treatment work records | Each format saved 8 of 8, comparator exit 0 | 58e220e8e03d7e4ce94d6f7557906f59fcd85c726f220b9dc93802ef7db6f252 |
| Plain untimed run on the exact timing inputs | All 11 Result fields equal for all 24 rows, 2,304 bytes per role | dfb8b0d490fee4c5cdcda1e3f4c6e51c452efa853b51f6d68c129a4f1f2f3ba9 |

The separate existing test_core and test_format3 fixtures both passed against the treatment header. The plain input FASTA has four contigs. Its 11-site reference and default candidate cap of 10,000 were identical in the plain parity and timing modes. The instrumented suite uses an extended 16-site reference and additional limits.

## Reproducing the comparison

From the repository root at the treatment revision, use a fresh temporary directory. The measured compiler was Darwin Command Line Tools clang++ with the MacOSX27.2 SDK and C++17/O2. The following shell commands reproduce the measurement layout. Each role uses one plain core header. Staging the driver and fixture together avoids resolving an archived sibling header. Check the emitted dependency files for the expected role paths. The original measurements ran serially with LC_ALL=C.

```sh
set -eu
export LC_ALL=C
work=$(mktemp -d)
mkdir "$work/driver" "$work/control" "$work/treatment"
git show 737ad24784ee9a2445a3eaf4f8262c87561a72f7:src/preservation_core.h > "$work/control/preservation_core.h"
cp src/preservation_core.h "$work/treatment/"
cp tests/native/bench_normalize_reuse.cpp tests/native/normalize_reuse_fixture.h "$work/driver/"
compiler=/Library/Developer/CommandLineTools/usr/bin/clang++
sdk=/Library/Developer/CommandLineTools/SDKs/MacOSX27.2.sdk
for role in control treatment; do
    "$compiler" -std=c++17 -O2 -isysroot "$sdk" -I "$work/$role" -MMD -MF "$work/$role.d" \
        "$work/driver/bench_normalize_reuse.cpp" -o "$work/$role-bench"
    "$work/$role-bench" --parity "$work/$role-parity-fixture" > "$work/$role.parity"
done
cmp "$work/control.parity" "$work/treatment.parity"
for run in c1 t1 t2 c2 c3 t3; do
    case "$run" in c*) role=control ;; t*) role=treatment ;; esac
    /usr/bin/time -l "$work/$role-bench" "$work/$run" > "$work/$run.tsv" 2> "$work/$run.time"
done
```

For the work diagnostic, use separate instrumented copies of both pinned headers and the same compiler/SDK flags. The hooks count entries to the actual normalize function. They are absent from the plain binaries.

```sh
cp tests/native/test_normalize_reuse.cpp "$work/driver/"
for role in control treatment; do
    mkdir "$work/$role-diag"
    awk -f tests/native/normalize_reuse_diag.awk "$work/$role/preservation_core.h" > "$work/$role-diag/preservation_core.h"
    "$compiler" -std=c++17 -O2 -isysroot "$sdk" -I "$work/$role-diag" -MMD -MF "$work/$role-diag.d" \
        "$work/driver/test_normalize_reuse.cpp" -o "$work/$role-diag/test"
    "$work/$role-diag/test" "$work/$role-diag/fixture" > "$work/$role-diag.parity"
done
cmp "$work/control-diag.parity" "$work/treatment-diag.parity"
"$work/treatment-diag/test" --compare-work "$work/control-diag/fixture/work.txt" "$work/treatment-diag/fixture/work.txt"
```

## Source and artifact identity

| Input or measured binary | SHA-256 |
| --- | --- |
| Control src/preservation_core.h | e165157b2fd046836b0b650940da5371296618f88d37f99250324d4785135300 |
| Treatment src/preservation_core.h | 42cbd86be7e8e85dd67b11713db3f3483997be96f6ace69b54d655d0ca7decf5 |
| bench_normalize_reuse.cpp | 05e4746b2e4a0db5d4045335f49830b4d6b1854387309a4e25c2b131fb4084ca |
| normalize_reuse_fixture.h | 4c62c19a42ea39efc63d7dbb96d8dd711996e0d1aaaae22b112d1b9ef10ea891 |
| normalize_reuse_diag.awk | 3fa306e95535071bef77e872a45518f8842f4fef983cdb2d76c0639e272e39e0 |
| test_normalize_reuse.cpp | 199aa4c117968a462de12b4ad915a5a1b040a1d6dbd0894f1c8fc47aa5403d61 |
| Compiler executable | fc5e3cbba6ad54059295cad716aacb05ace8235f5f7a6031f1d8c15d24bd6d66 |
| SDKSettings.json | 82165b15ee21419d4cffaba8895797391b908ac9499ba70f561f3e1bd16758d6 |
| Measured control binary | f3613882061b2726aebd95c312c95023cbcd23e0a4915caf571e62a8059d7d89 |
| Measured treatment binary | 0c87503aa5433e2c185ce89eb9a5af760085188c6cab9a205eb4634dec5e54f9 |

SHA-256 checksums of the archived raw measurements follow in execution order.

```text
b30cddd4228266196764782ded3ed66fb53ad71f7fddcf71b067d8dc4f7749e6  c1.tsv
cb299f65173a1c30ad2dbaa595077429b7c674e398361989197614322a6b703f  c1.time
f1524c241a2e8de0d784ac69831c7acf0607ff07c38fd75e776e0a5b306e4749  t1.tsv
df3a9d4242554a8089e97491dde1edab671e916a08b34f26fda5b9e4bb705205  t1.time
6b995b2068cdaedcac49fcb7e957a5fb26af2ef10fae0458a501cb06e864705d  t2.tsv
cb24610f09782e5fdceaee94bf6ef30f3554c4a859047505e80f6af2d56bec20  t2.time
e6342facaea11dae20a1fc709c762acc4de03f4f45cd849d165e4e1e8a86c814  c2.tsv
89556423e0afde3fa56b4dc790b910b5e0637a6ec8cdb6a10edb23ce83aec259  c2.time
003e84b80cc918f5bdea97289a625cfc640addc17e7d15e896ca4fc6aa0b3915  c3.tsv
d5fb48fdc6ed1f4d254a5de2a58b30cb982858e771b45dd75831d35cf325734c  c3.time
b5142541c779e315af38c2125e3acf5e7884435a62f373e25ad0af8c96a98e87  t3.tsv
02f30b20c5abde59847882713108360bc4e7d3c8a41ddc072d5c2b47ce73219a  t3.time
```
