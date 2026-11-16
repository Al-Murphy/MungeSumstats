"""Retain every sanitizer stderr and fail on any diagnostic, using built binaries."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
OUT = HERE / 'sanitized-strict'
OUT.mkdir()
os.environ['ASAN_OPTIONS'] = 'halt_on_error=1'
os.environ['UBSAN_OPTIONS'] = 'halt_on_error=1:print_stacktrace=1'
spec = importlib.util.spec_from_file_location('unchanged_count_tests', HERE / 'test_count.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.HERE = HERE / 'sanitized'
original_run = module.run
calls = []


def checked_run(binary, *args):
    result = original_run(binary, *args)
    calls.append({'binary': binary, 'stderr_hex': result.stderr.hex(), 'exit_code': result.returncode})
    assert not result.stderr, 'Unexpected sanitizer/writer diagnostic'
    return result


module.run = checked_run
count = 0
for version in [2, 3]:
    for case, rows in module.CASES.items():
        with tempfile.TemporaryDirectory(dir=OUT) as temp:
            module.test_exact_native_sizes(Path(temp), version, case, rows)
            count += 1
with tempfile.TemporaryDirectory(dir=OUT) as temp:
    module.test_format_suffix_refused(Path(temp))
    count += 1
arithmetic = subprocess.run([str(module.HERE / 'arithmetic_fixture')], capture_output=True)
assert arithmetic.returncode == 0 and not arithmetic.stderr
assert len(calls) == 2 * (count - 1) and all(not row['stderr_hex'] for row in calls)
(OUT / 'summary.json').write_text(json.dumps({
    'status': 'success', 'unchanged_differential_cases': count, 'native_calls': calls,
    'arithmetic_fixture': {'exit_code': arithmetic.returncode, 'stderr_hex': arithmetic.stderr.hex()},
    'sanitizers': ['address', 'undefined'], 'halt_on_error': True,
    'source_test_sha256': hashlib.sha256((HERE / 'test_count.py').read_bytes()).hexdigest(),
    'earlier_run_limitation': 'Earlier runner checked exit statuses but did not inspect successful subprocess stderr; this fresh strict run closes that harness gap',
}, indent=2) + '\n')
print(json.dumps({'status': 'success', 'cases': count, 'successful_native_calls_with_empty_stderr': len(calls)}))
