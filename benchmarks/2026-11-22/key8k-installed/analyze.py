"""Independently reconcile returned metadata, then summarize matched pairs."""
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

OP = Path('/private/tmp/mss-combined-operator-20261122')
HERE = Path(__file__).resolve().parent
RET = OP / 'returned'

def read(path):
    return json.loads(path.read_bytes())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def tsv(path):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    return rows

assert read(OP / 'terminal.json')['exit_code'] == 0
assert read(OP / 'copy-terminal.json')['exit_code'] == 0
copy = read(OP / 'controller-copy-proof.json')
assert not copy['absent'] and len(copy['controller_pins']) == 88
assert sum(pin['bytes'] for pin in copy['controller_pins'].values()) == 905074
for name, pin in copy['controller_pins'].items():
    path = RET / name
    assert path.stat().st_size == pin['bytes'] and sha(path) == pin['sha256'], name
pins = read(OP / 'PAYLOAD-PINS.json')
assert {name: pin['sha256'] for name, pin in pins.items()} == copy['actual_source_pins']
proof = read(RET / 'control/combined-proof.json')
assert sha(RET / 'control/combined-proof.json') == '69b58946bacd0e4160bfde0964368b486c99216209fa19f9306b5ad905d01a07'
assert proof['status'] == 'combined_installed_R_guard_complete'
assert proof['no_signaled_or_unknown_cleanup'] and proof['all53_retained_identity_sha_unchanged']
assert proof['stable_file_comparisons'] == 112 and proof['fresh_workers'] == 12
assert proof['independent19_field_oracle_rows'] == [1000000, 10000000]
assert proof['final_storage']['peak_setup_growth_bytes'] < 80000000000
assert proof['final_storage']['workspace_conservative_free_bytes'] > 186 * 1024**3
receipt = read(RET / 'output/comparison-receipt.json')
assert receipt['status'] == 'combined_installed_R_scientifically_verified'
assert sha(RET / 'output/comparison-receipt.json') == '86f1f63121ac552ccbfaae62b5b01317592c5178332e047f33b3fced898d094f'
assert receipt['recipe']['candidate_commit'] == proof['candidate_commit'] == '8bc67611e96c588db79e853271951e2611056756'
qualified = read(RET / 'output/functional/verified.json')
cases = tsv(RET / 'output/functional/receipt/tests.tsv')
assert len(cases) == qualified['tests'] == 269
assert sum(int(row['passed']) for row in cases) == qualified['passed'] == 3077
assert all(row['skipped'] == row['error'] == 'FALSE' and int(row['failed']) == 0 for row in cases)
assert sum(int(row['warning']) for row in cases) == qualified['test_warnings'] == 2
assert tsv(RET / 'output/functional/receipt/conditions.tsv') == qualified['conditions']
assert [(r['file'], r['test'], int(r['passed'])) for r in cases] == [
    (r['file'], r['test'], r['passed']) for r in qualified['cases']]
assert receipt['qualified_candidate'] == qualified

def number(value):
    result = float(value)
    assert math.isfinite(result)
    return result

def gc_cpu(metrics, phase):
    return sum(number(metrics[f'{phase}.{i}']) for i in (1, 2, 4, 5))

fields = {
    'operation_seconds': lambda m, p: number(m['processing.elapsed']),
    'operation_plus_closing_gc_seconds': lambda m, p: number(m['processing_plus_closing_gc.elapsed']),
    'closing_gc_seconds': lambda m, p: number(m['closing_gc.elapsed']),
    'worker_wall_seconds': lambda m, p: number(p['wall_seconds']),
    'worker_cpu_seconds': lambda m, p: number(p['user_seconds']) + number(p['system_seconds']),
    'peak_rss_bytes': lambda m, p: p['peak_rss_bytes'],
    'operation_gc_cpu_seconds': lambda m, p: gc_cpu(m, 'gc_operation_delta'),
    'operation_gc_elapsed_seconds': lambda m, p: number(m['gc_operation_delta.3']),
    'operation_plus_closing_gc_cpu_seconds': lambda m, p: gc_cpu(m, 'gc_operation_plus_closing_delta'),
    'operation_plus_closing_gc_collection_elapsed_seconds': lambda m, p: number(m['gc_operation_plus_closing_delta.3']),
    'gc_through_pre_receipt_cpu_seconds': lambda m, p: gc_cpu(m, 'gc_through_receipt_delta'),
    'gc_through_pre_receipt_elapsed_seconds': lambda m, p: number(m['gc_through_receipt_delta.3']),
}
primary = tuple(fields)
for phase in ('gc_setup_delta', 'gc_reference_delta', 'gc_preparation_delta',
              'gc_operation_delta', 'gc_closing_delta', 'gc_operation_plus_closing_delta',
              'gc_through_receipt_delta'):
    for index, name in enumerate(('user_cpu', 'system_cpu', 'elapsed', 'child_user_cpu', 'child_system_cpu'), 1):
        fields[f'{phase}.{name}'] = lambda m, p, ph=phase, i=index: number(m[f'{ph}.{i}'])

summary = {'candidate': proof['candidate_commit'], 'baseline': 'qualified R-lifetime-only package; 16KiB key pages',
           'candidate_change': '8KiB key pages at the same 8MiB payload budget',
           'negative_percent_is_improvement': True, 'results': {}}
for label, rows in (('1M', 1000000), ('10M', 10000000)):
    scale = receipt['results'][label]
    assert scale['paired_order'] == ['AB', 'BA', 'AB'] and scale['rows'] == rows
    assert scale == read(RET / f'output/verification/{label}/scale-receipt.json')
    stable = read(RET / f'output/verification/{label}/stable-editions.json')
    digests = read(RET / f'output/verification/{label}/artifact-digests.json')
    assert len(digests) == stable['edition_count'] == 7 and stable['stable_byte_equality']
    assert len(stable['normalized_columns']) == 19 and stable['normalized_rows_checked'] == rows
    maps = list(digests.values())
    names = set(maps[0]) - {'manifest.rds'}
    assert len(names) == 8
    for name in names:
        assert all(d[name] == maps[0][name] for d in maps)
    assert (RET / f'output/verification/{label}/stable-manifests.txt').read_bytes() == b'pass\n'
    oracle = read(RET / f'output/verification/{label}/whole-SNP-oracle.json')
    assert oracle['passed'] and oracle['all_rows_examined'] and oracle['unaccounted_deltas'] == 0
    assert oracle['rows_expected'] == oracle['versions']['mss7']['counts']['all_19_fields_correct_rows'] == rows
    records = {}
    namespaces = []
    for record in scale['measurements']:
        name = record['run']
        base = RET / f'output/runs/{label}/{name}'
        metrics = {r['key']: r['value'] for r in tsv(base / 'metrics.tsv')}
        process = read(base / 'process.json')
        assert record['metrics'] == metrics and record['process'] == process
        assert process['exit_code'] == 0 and metrics['threads'] == '1'
        assert number(metrics['source_rows']) == number(metrics['normalized_rows']) == rows
        assert record['role'] == metrics['role'] == name[-1]
        roster = tsv(base / 'namespaces.tsv')
        namespaces.append([{k: v for k, v in r.items() if k != 'path'} if r['package'] == 'MungeSumstats' else r for r in roster])
        records[name] = (metrics, process)
    assert len(records) == 6 and all(n == namespaces[0] for n in namespaces)
    comparisons = []
    for pair in (1, 2, 3):
        a, b = records[f'pair-{pair:02d}-A'], records[f'pair-{pair:02d}-B']
        values = {}
        for name, extract in fields.items():
            av, bv = extract(*a), extract(*b)
            values[name] = {'baseline': av, 'candidate': bv, 'difference': bv-av,
                            'percent': 100*(bv-av)/av if av else None}
        comparisons.append({'pair': pair, 'order': ('AB', 'BA', 'AB')[pair-1], 'values': values})
    medians = {}
    for name in fields:
        values = [c['values'][name] for c in comparisons]
        percents = [v['percent'] for v in values if v['percent'] is not None]
        medians[name] = {'baseline': statistics.median(v['baseline'] for v in values),
                         'candidate': statistics.median(v['candidate'] for v in values),
                         'paired_difference': statistics.median(v['difference'] for v in values),
                         'paired_percent': statistics.median(percents) if percents else None}
    summary['results'][label] = {'pairs': comparisons, 'medians': medians}

with (HERE / 'analysis.json').open('x') as stream:
    json.dump(summary, stream, indent=2, allow_nan=False)
    stream.write('\n')
print(json.dumps({label: {name: value for name, value in r['medians'].items() if name in primary}
                  for label, r in summary['results'].items()}, indent=2))
