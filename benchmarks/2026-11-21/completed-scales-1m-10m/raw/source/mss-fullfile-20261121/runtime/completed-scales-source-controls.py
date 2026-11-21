"""Source-only controls for the scoped verifier."""
from copy import deepcopy
import json
from pathlib import Path
import runpy
from types import SimpleNamespace
from unittest.mock import patch


source = Path('/private/tmp/mss-fullfile-20261121/runtime/verify-completed-scales.py')
compile(source.read_bytes(), str(source), 'exec')
print('syntax: pass')
module = runpy.run_path(str(source))
warm = module['load_frozen']()
print('authenticated source imports: pass')


def rejects(expected, function, *args):
    try:
        function(*args)
    except ValueError as exc:
        assert expected in str(exc), str(exc)
    else:
        raise AssertionError(f'Expected rejection: {expected}')


controls = warm.Controls()
config = module['read_pinned'](controls, module['MATRIX'] / 'pins.json', module['OPENING_SHA256'])['configuration']
roster = warm.d.expected_roster()
cases = module['select_cases'](config, roster)
assert [(size, case['rows'], case['pairs']) for size, case in cases.items()] == [('1M', 1_000_000, 3), ('10M', 10_000_000, 3)]
incomplete = deepcopy(config)
incomplete['cases'][0]['pairs'] = 2
rejects('Unexpected selected case: matched-1M', module['select_cases'], incomplete, roster)
print('in-memory roster selection and wrong pair count rejection: pass')

failure, _ = module['failure_evidence'](warm, controls)
assert failure['returned_outcome']['raw_exit'] == -15
bad_failure = {**failure, 'terminal_receipt_present': True}
failed_bytes = json.dumps(bad_failure).encode()
with patch.object(controls, 'data', side_effect=lambda path: (
        failed_bytes if path == module['FAILURE'] else warm.d.control_bytes(path))):
    rejects('Original guard failure binding differs', module['failure_evidence'], warm, controls)
print('failure binding and false terminal record rejection: pass')

with patch.dict(module['worker_dirs'].__globals__, {'MATRIX': source / 'synthetic-missing-directory'}):
    rejects('Missing or indirect retained directory', module['worker_dirs'], cases['1M'], warm.r.ROLES)
print('missing worker rejection: pass')
rejects('Unexpected control digest:', module['read_pinned'], controls, module['CONFIG'], '0' * 64)
print('wrong configuration digest rejection: pass')

# Exercise the real collector with one in-memory worker and record requested pins.
dest = Path('/synthetic/pair-01-baseline')
process = {'started_utc': '2026-11-21T08:14:10+00:00', 'finished_utc': '2026-11-21T08:14:11+00:00'}
metrics = {f'{phase}.{stat}': '1' for phase in ('package_setup', 'reference_open', 'processing')
           for stat in ('elapsed', 'user.self', 'sys.self', 'user.child', 'sys.child')}
requested = []
inventory = {'assets': {}, 'dependencies': {}, 'r_toolchain': {}}
with patch.object(warm, 'checked_run', return_value=({}, process, metrics, {})):
    module['collect_workers'](warm, SimpleNamespace(pin=requested.append), config, {'1M': cases['1M']},
                              {'1M': [(1, 'baseline', dest)]}, inventory,
                              {'started_utc': '2026-11-21T08:14:09+00:00'}, failure)
missing = {dest / name for name in ('process.stdout.log', 'process.stderr.log', 'r-heap.rds')} - set(requested)
assert not missing, f'Missing worker log pins: {sorted(map(str, missing))}'
print('worker log pins include the heap diagnostic: pass')
