"""Run the unchanged29 differential assertions against separate sanitized binaries."""
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
CORE = Path('/private/tmp/gwas-harvest-stage2/mss-context-format3-release/src/preservation_core.h')
OUT = HERE / 'sanitized'
OUT.mkdir()
pins = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [CORE, HERE / 'count_reference.cpp', HERE / 'writer_fixture.cpp',
         HERE / 'arithmetic_fixture.cpp', HERE / 'test_count.py']}
for name in ['count_reference', 'writer_fixture', 'arithmetic_fixture']:
    result = subprocess.run(['c++', '-std=c++17', '-O1', '-g',
                             '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                             '-I', str(CORE.parent), str(HERE / (name + '.cpp')),
                             '-o', str(OUT / name)], capture_output=True, text=True)
    (OUT / (name + '-build.log')).write_text(result.stdout + result.stderr)
    result.check_returncode()
spec = importlib.util.spec_from_file_location('unchanged_count_tests', HERE / 'test_count.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.HERE = OUT  # Only binary location differs; every original assertion is unchanged.
count = 0
for version in [2, 3]:
    for case, rows in module.CASES.items():
        with tempfile.TemporaryDirectory(dir=OUT) as temp:
            module.test_exact_native_sizes(Path(temp), version, case, rows)
            count += 1
with tempfile.TemporaryDirectory(dir=OUT) as temp:
    module.test_format_suffix_refused(Path(temp))
    count += 1
arithmetic = subprocess.run([str(OUT / 'arithmetic_fixture')], capture_output=True, text=True)
assert arithmetic.returncode == 0 and not arithmetic.stderr
for name, expected in pins.items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == expected
(OUT / 'summary.json').write_text(json.dumps({
    'status': 'success', 'unchanged_differential_cases': count,
    'arithmetic_fixture': 'passed', 'sanitizers': ['address', 'undefined'],
    'scope': 'macOS/Apple compiler; no nativeWindows claim; no full/reference build',
    'source_pins': pins,
}, indent=2) + '\n')
print(json.dumps({'status': 'success', 'cases': count, 'arithmetic_fixture': 'passed'}))
