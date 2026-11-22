"""Independent arithmetic over unchanged, returned native benchmark metadata."""
from pathlib import Path
from decimal import Decimal
from statistics import median
import hashlib
import json

HERE = Path(__file__).resolve().parent
OP = Path('/private/tmp/mss-batch-rank-benchmark-operator-fix1-20261122')
ROOT = OP / 'returned'
GOLD = Path('/private/tmp/mss-key-page-8k-experiment-20261122/expected')
pins = json.loads((OP / 'returned-pins.json').read_bytes())
assert len(pins) == 92
for name, pin in pins.items():
    raw = (ROOT / name).read_bytes()
    assert len(raw) == pin['bytes']
    assert hashlib.sha256(raw).hexdigest() == pin['sha256']
assert (ROOT / 'source-before.txt').read_bytes() == (ROOT / 'source-after.txt').read_bytes()
assert all(line.endswith(': OK') for line in (ROOT / 'source-before.txt').read_text().splitlines())

stored = json.loads((ROOT / 'paired.json').read_text(), parse_float=Decimal)
schedule = json.loads((GOLD / 'order.json').read_bytes())
out = {}
for order in ('matched', 'shuffled'):
    pairs = []
    for number in (1, 2, 3):
        values = {}
        for role in ('baseline', 'candidate'):
            stem = ROOT / f'panel-{order}-{number}-{role}'
            assert stem.with_suffix('.json').read_bytes() == (GOLD / 'panel.json').read_bytes()
            t = json.loads(stem.with_suffix('.timing.json').read_bytes())
            ns = t.pop('harmonise_ns')
            batch = t.pop('batch_rank', None)
            assert t == schedule['panel'][order]
            if role == 'candidate':
                assert batch['attempts'] == batch['admitted'] == 10
                assert batch['arena_peak_bytes'] <= 32 * 1024**2
                assert batch['gather_max_keys'] <= 200000
                assert batch['gather_max_key_bytes'] <= 224
            else:
                assert batch is None
            p = json.loads(stem.with_suffix('.process.json').read_text(), parse_float=Decimal)
            assert p['exit_code'] == 0 and p['pid'] > 0
            assert p['command'] == [f'/workspace/mss-batch-rank-benchmark-results-20261122/{role}',
                                    'panel', '/workspace/mss-task4/inputs/reference/matched-index',
                                    '/workspace/mss-task4/inputs/reference/GRCh37.ensembl104.chr22.fa', order]
            values[role] = dict(native_seconds=Decimal(ns) / 10**9,
                                wall_seconds=p['wall_seconds'],
                                cpu_seconds=p['user_seconds'] + p['system_seconds'],
                                peak_rss_bytes=Decimal(p['peak_rss_bytes']))
            assert p['wall_seconds'] >= values[role]['native_seconds'] > 0
        delta = {k: values['candidate'][k] - values['baseline'][k] for k in values['baseline']}
        percent = {k: 100 * delta[k] / values['baseline'][k] for k in delta}
        original = stored[order]['pairs'][number - 1]
        assert original['execution_order'] == (['candidate', 'baseline'] if number == 2 else ['baseline', 'candidate'])
        for field in percent:
            assert abs(percent[field] - original['signed_percent'][field]) < Decimal('0.00000000001')
        assert percent['native_seconds'] > 0
        pairs.append(dict(pair=number, metrics=values, signed_delta=delta, signed_percent=percent))
    medians = {field: median(p['signed_percent'][field] for p in pairs) for field in pairs[0]['signed_percent']}
    for field, value in medians.items():
        assert abs(value - stored[order]['median_signed_percent'][field]) < Decimal('0.00000000001')
    out[order] = dict(pairs=pairs, median_signed_percent=medians)

closure_root = Path('/private/tmp/mss-batch-rank-benchmark-closure-20261122')
closure_raw = (closure_root / 'closure.stdout').read_bytes()
assert hashlib.sha256(closure_raw).hexdigest() == 'dc7beb39b434ded686030652cc98ad994293c751fa43993fe3aa32c969bb884e'
closure = json.loads(closure_raw)
assert json.loads((closure_root / 'terminal.json').read_bytes())['exit_code'] == 0
assert closure['actual_native_child_rc'] == 0 and closure['actual_unsignaled_quiescent_drain'] is True
assert closure['benchmark_rerun'] is False and closure['metadata_pins'] == pins
assert closure['compiler_at_closing']['sha256'] == '1353e9bdd29a7295c7226bf6c63abccce056d8cac31f112e5cdbecc3f28c2769'
assert closure['paired'] == json.loads((ROOT / 'paired.json').read_bytes())
record = dict(scope='native repeated-SNP 1M panel, three alternating pairs per query order',
              decision='Do not integrate the batch-rank prototype into the R API',
              source_manifest_sha256='e87a5cd8883074529c4146711bb353dfc64008d2a9013a923b50a7c5708e1c7e',
              candidate_header_sha256='9f71fa92eaa8e6b8d249ba58d7863630a7f6958e491ea413a48d521151932225',
              baseline_header_sha256='5881233fbad917d4475f4801497c90a720428c41ba5afb8a7d38c234da60f155',
              returned_metadata_files_verified=92, complete_panel_digests_verified=12,
              all_six_native_pairs_slower=True, original_guard_exit_code=1,
              original_guard_failure='tuple/list mismatch in final receipt comparison',
              verification_only_completion_recorded_separately=True,
              verification_only_proof_sha256=hashlib.sha256(closure_raw).hexdigest(),
              unrecoverable_guard_local_fields=closure['unrecoverable_volatile_preclosing_fields'],
              results=out)
target = HERE / 'analysis.json'
target.write_text(json.dumps(record, indent=2, default=str) + '\n')
print(json.dumps({k: {m: str(v) for m, v in x['median_signed_percent'].items()} for k, x in out.items()}, indent=2))
