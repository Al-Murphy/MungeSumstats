import csv
import gzip
import hashlib
import io
import json
import math
import re
import statistics
from collections import Counter
from pathlib import Path

T = Path('/private/tmp')
R = T/'mss-fullfile-20261121/runtime'
B = T/'mss-key-geometry-builds-20261121'
D = T/'mss-key-geometry-benchmark-20261121/diagnostic'
O = T/'mss-key-geometry-oracle-20261121'
P = T/'mss-key-geometry-preparation-20261121/metadata'
H = T/'mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade'
baseline = 'f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b'
candidate = '8b85aa5ffc45ad7632cfc5319e85dd5eb3089be6'
commits = {'baseline': baseline, 'candidate': candidate}
expected = {'hard-cases-f2': (32, 2), 'hard-cases-f3': (32, 3),
            'declared-f2': (9, 2), 'declared-f3': (9, 3), 'matched-1M': (1000000, 2)}
order = ['pair-01-baseline', 'pair-01-candidate', 'pair-02-candidate',
         'pair-02-baseline', 'pair-03-baseline', 'pair-03-candidate']
fields = set('source_row status chr pos ref alt id effect eaf swapped complemented candidate_count entity_type entity_id normalized_build grch38_chr grch38_pos grch38_ref grch38_alt'.split())
phases = ('package_setup', 'reference_open', 'processing')
phase_components = ('elapsed', 'user.self', 'sys.self', 'user.child', 'sys.child')

def read(path):
    with Path(path).open('rb') as stream:
        payload = stream.read(600001)
    assert len(payload) <= 600000, path
    return payload

def digest(data):
    return {'bytes': len(data), 'md5': hashlib.md5(data).hexdigest(), 'sha256': hashlib.sha256(data).hexdigest()}

def sha(path):
    return digest(read(path))['sha256']

def j(path):
    return json.loads(read(path))

def pin(record):
    assert sha(record['path']) == record['sha256'], record
    return Path(record['path'])

def table(path, header):
    rows = list(csv.reader(io.StringIO(read(path).decode()), delimiter='\t', strict=True))
    assert rows.pop(0) == header
    assert all(len(row) == len(header) for row in rows)
    keys = {row[0]: row[1:] for row in rows}
    assert len(keys) == len(rows)
    return keys

def kv(path, header=('key', 'value')):
    return {k: v[0] for k, v in table(path, list(header)).items()}

def normalized(path):
    incoming = csv.DictReader(io.StringIO(read(path).decode()), delimiter='\t', strict=True)
    assert len(incoming.fieldnames) == len(fields) and set(incoming.fieldnames) == fields
    rows = list(incoming)
    assert all(set(row) == fields and None not in row.values() and row['source_row'] == str(i)
               for i, row in enumerate(rows, 1))
    return rows

def guard(record, argv=None):
    planned, running, terminal = [json.loads(line) for line in read(pin(record)).splitlines()]
    assert planned['status'] == 'planned' and planned['allowed_exits'] == [0] and running['status'] == 'running'
    assert planned['label'] == running['label'] == terminal['label']
    assert planned['utc'] <= running['utc'] <= terminal['utc'] and planned['receipt'] == record['path']
    assert running['pid'] == running['expected_pgid'] == running['expected_sid'] > 1
    assert terminal['status'] == 'completed' and terminal['raw_exit'] == 0 and terminal['closure'] == 'confirmed'
    assert terminal['leader_reaped'] and not terminal['cleanup_errors'] and not terminal['signals']
    assert argv is None or planned['argv'] == argv
    return planned, running, terminal

def observed(label, expected_hash, outer=True):
    path = T/f'mss-key-geometry-{label}-observation-20261121.json'
    assert sha(path) == expected_hash
    result = j(path)
    assert not result['inspection_errors'] and not result['live_members']
    if outer:
        assert result['actual_exit_code'] == 0 and result['actual_terminal_tool_chunk']
    else:
        assert 'actual_exit_code' not in result and 'actual_terminal_tool_chunk' not in result
    return result

def attested_guard(result, argv=None):
    planned, running, terminal = guard(result['guard'], argv)
    assert result['terminal'] == terminal and result['absent_group_session'] == running['pid']
    return planned, running, terminal

bo = observed('builds', 'ee38ae6c506f2a700e3159a0d6ce038d3785a32f634e00dee7de72b03ec39473')
op = j(pin(bo['operator_result']))
assert bo['producer'] == op['producer']
pin(bo['producer'])
pin(op['reused_producer'])
assert op['status'] == 'completed' and op['commits'] == op['validated_commits'] == bo['commits'] == commits
assert len(op['phases']) == len(bo['phases']) == 9
builds = {}
for phase, observation in zip(op['phases'], bo['phases']):
    assert observation['guard'] == phase['receipt']
    plan, run, end = guard(observation['guard'], phase['argv'])
    assert plan['label'] == phase['label'] == observation['label'] and end == phase['terminal']
    assert run['pid'] == observation['group_session'] and observation['group_session'] in bo['absent_groups_sessions']
    assert observation['closure'] == 'confirmed' and observation['leader_reaped'] and observation['raw_exit'] == 0
    for p in phase['transcripts']:
        if Path(p['path']).suffix != '.tar':
            pin(p)
provenance = j(pin(op['provenance']))
for role, commit in commits.items():
    path = pin(op['roles'][role])
    build = builds[role] = j(path)
    assert provenance[role] == op['roles'][role] and build['source_commit'] == commit
    library = B/role/'library'
    assert build['library'] == str(library) and {x.name for x in library.iterdir()} == {'MungeSumstats'}
    assert build['dll']['path'] == str(library/'MungeSumstats/libs/MungeSumstats.so')
    assert build['dll_sha256'] == build['dll']['sha256']
    assert read(B/role/'source-tree.stdout').decode().strip() == build['source_tree']
    assert build['install_argv'][-2:] == [f'--library={library}', build['archive']['path']]
    planned, running, terminal = guard(build['install_transcripts'][0], build['install_argv'])
    assert terminal == build['terminal']
    stderr = read(pin(build['install_transcripts'][2])).decode()
    assert all(s in stderr for s in ('** using staged installation', 'using C++17', 'temporary location',
                                     'final location', '* DONE (MungeSumstats)'))
    print('BUILD', role, build['source_tree'], build['dll_sha256'], planned['label'])
assert len({build['source_tree'] for build in builds.values()}) == 2
assert len({build['dll_sha256'] for build in builds.values()}) == 2
assert builds['baseline']['compile_and_link_argv'] == builds['candidate']['compile_and_link_argv']
assert len(builds['baseline']['compile_and_link_argv']) == 3
assert builds['baseline']['toolchain'] == builds['candidate']['toolchain']

admission = j(D.parent/'admission.json')
for name in ('build_observation', 'preparation_observation', 'preparation_binding', 'configuration'):
    pin(admission[name])
assert (admission['baseline_commit'], admission['candidate_commit']) == (baseline, candidate)
assert not admission['competing_processes'] and admission['available_bytes'] >= admission['disk_bytes'] + admission['floor_bytes']
assert admission['disk_bytes'] == 3*1024**3 and admission['floor_bytes'] == 186*1024**3
bmo = observed('benchmark', '004e670d742c268adf2ca342b90cd41b6e34c85e338de47796dba6b0a93b5775')
bmplan, _, bmterm = attested_guard(bmo, admission['argv'])
assert bmterm['minimum_free_bytes'] >= admission['floor_bytes'] and not bmterm['missed_samples']
assert bmterm['peak_disk_allocated_bytes'] < admission['disk_bytes']
assert bmterm['peak_sampled_group_rss_bytes'] < admission['rss_bytes']
assert bmterm['elapsed_seconds'] < admission['wall_seconds']
assert not read(bmplan['stdout']) and not read(bmplan['stderr'])
po = observed('preparation', '04d185cc80ac89f004c90ea2fc1afc8f809d21898e09a4044894c313b64fe491')
pp, _, pe = attested_guard(po)
assert pp['argv'][1:6] == ['-I', '-B', '-S', str(R/'key-geometry-bench.py'), 'prepare']
assert pp['argv'][-2:] == ['--candidate', candidate]
assert pe['utc'] < bmplan['utc']

summary = j(pin(bmo['summary']))
policy = j(pin(bmo['namespace_policy']))
binding = j(pin(policy['preparation']))
cfg = j(pin(policy['configuration']))
assert admission['preparation_binding'] == policy['preparation']
assert admission['configuration'] == policy['configuration'] == binding['configuration']
assert binding['provenance'] == op['provenance'] and binding['candidate'] == candidate
assert policy['source_commits'] == commits and not policy['required_baseline_only']
for name in ('binding', 'bootstrap', 'baseline_override', 'frozen_sources', 'matrix_config', 'key_geometry_override'):
    assert policy[name] == binding[name]
for name in ('binding', 'bootstrap', 'matrix_config'):
    pin(policy[name])
for p, control in policy['frozen_sources'].items():
    assert sha(p) == control
assert policy['binding']['sha256'] == 'a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94'
over = policy['key_geometry_override']
pin(over['wrapper'])
pin(over['counter_review'])
assert over['wrapper']['sha256'] == admission['wrapper_sha256'] == 'c2eb2af6b254a981272740eb05b6632b4dc0b55603f62c6cb9a1149a62dee976'
assert over['counter_review']['sha256'] == 'd478ab3a248afb1f94cca78c913037735aadea29d8e6bd9c4fbbd539c4f425d0'
assert over['actual_baseline'] == policy['baseline_override']['actual'] == baseline
assert over['original_baseline'] == '14ea7af51a19272c7792fbe582a3d42b291b4a1a'
assert len(over['prior_eager_only']) == 13 and over['revised_eager_only'] == []
assert list(summary['cases']) == list(policy['cases']) == list(expected)
assert summary['fresh_processes'] and summary['suite'] == 'diagnostic' and summary['profile'] == 'none'
assert summary['target_status'] == 'not_applicable' and summary['processing_speed_and_rss_target_met'] is None
assert bmplan['utc'] <= summary['started_utc'] < summary['finished_utc'] <= bmterm['utc']

inventory = read(pin(binding['preparation_inventory']))
assert inventory == read(D/'namespace-inputs-before.json') == read(D/'namespace-inputs-after.json')
iv = json.loads(inventory)
libraries = iv['installed_dependency_tree_sha256']
assert sorted(map(len, libraries.values())) == [30,96]
pre, post = j(D/'pins.json'), j(D/'post-pins.json')
assert pre['configuration_file'] == {'path':policy['configuration']['path'], **digest(read(policy['configuration']['path']))}
assert pre['assets'] == post['assets'] and pre['installed_packages'] == post['installed_packages'] == iv['role_trees']
for name, p in pre['harness'].items():
    assert digest(read(H/name)) == p
print('PINS', len(pre['assets']), [len(v) for v in libraries.values()], digest(inventory)['sha256'])

cases = {case['id']: case for case in pre['configuration']['cases']}
assert list(cases) == list(expected)
ns = {(row['case'], row['worker']): row for row in policy['workers']}
osw = {(row['case'], row['worker']): row for row in bmo['workers']}
roster = [(case, worker) for case in expected for worker in order]
assert list(ns) == list(osw) == roster and len(policy['workers']) == len(bmo['workers']) == len(roster) == 30
assert len({x['pid'] for x in osw.values()}) == 30
shared = {}
last = summary['started_utc']
tiny_cells = 0
for case_id, (n, format_version) in expected.items():
    tiny = n < 100
    case = cases[case_id]
    assert case['rows'] == n and case['pairs'] == 3 and case['target'] is False
    ref = pre['configuration']['references'][case['reference']]
    assert format_version == ref['format']
    assert cfg['references'][case['reference']] == ref
    s = summary['cases'][case_id]
    assert s['rows'] == n and s['target'] is False
    assert s['cardinality'] == pre['cardinalities'][case_id]
    assert len(s['observations']) == 6
    v = j(D/case_id/'verification.json')
    assert s['verification'] == v and v['edition_count'] == 6
    assert v['stable_byte_equality'] and v['stable_rds_manifest_equality']
    assert v['normalized_rows_checked'] == n and sum(v['status_counts'].values()) == n
    assert len(v['normalized_columns']) == 19 and set(v['normalized_columns']) == fields
    assert v['oracle_rows_checked'] == (n if tiny else 0)
    assert set(v['oracle_fields'] or []) == (fields if tiny else set())
    gold = normalized(pin(case['oracle'])) if tiny else None
    if tiny:
        assert len(gold) == n
        original = read(pin(case['input']))
        decoded = gzip.decompress(original) if original[:2] == b'\x1f\x8b' else original
        assert digest(original) == v['source_artifact'] and digest(decoded) == v['decoded_source']
    r_cmd = [cfg['r_executable']['path'], '--vanilla', str(H/'worker.R')]
    m = j(D/case_id/'manifest-verification.json')
    assert m['exit_code'] == 0 and m['rss_complete'] and m['rss_source'] == 'wait4 ru_maxrss'
    assert m['command'] == r_cmd + ['manifests',str(D/case_id/'stable-manifests.txt')] + [
        str(D/case_id/worker/'edition/manifest.rds') for worker in order]
    assert read(D/case_id/'stable-manifests.txt') == b'pass\n'
    assert not read(D/case_id/'manifest-verification.stderr.log')
    assert not (D/case_id/'byte-deltas.json').exists()
    cert = policy['cases'][case_id]
    assert cert['baseline_only'] == cert['candidate_only'] == []
    artifacts = j(D/case_id/'artifact-digests.json')
    assert list(artifacts) == [str(D/case_id/worker) for worker in order]
    first = next(iter(artifacts.values()))
    assert first['source-artifact'] == v['source_artifact'] == pre['assets'][case['input']['path']]
    assert first['source.tsv'] == v['decoded_source']
    roles = {role: [] for role in commits}
    masked = []
    for worker, entry in zip(order, s['observations']):
        work = D/case_id/worker
        edition = work/'edition'
        role = worker.split('-')[-1]
        inp = kv(work/'settings.tsv')
        metrics = kv(work/'metrics.tsv')
        process = j(work/'process.json')
        assert entry['directory'] == str(work) and entry['pair'] == int(worker[5:7])
        assert entry['process'] == process and entry['metrics'] == metrics
        assert process['case'] == case_id and process['role'] == role
        assert process['pid'] == osw[case_id, worker]['pid']
        assert process['exit_code'] == osw[case_id, worker]['exit_code'] == 0
        assert process['rss_complete'] and process['rss_source'] == 'wait4 ru_maxrss' and process['peak_rss_bytes'] > 0
        assert process['command'] == r_cmd + ['run', str(work/'settings.tsv')]
        assert last <= process['started_utc'] < process['finished_utc']
        last = process['finished_utc']
        assert not read(work/'process.stderr.log')
        assert metrics['threads'] == '1' and float(metrics['source_rows']) == float(metrics['normalized_rows']) == n
        assert metrics['implementation'] == 'mss-preservation-7' and 'darwin' in metrics['R_platform']
        assert all(math.isfinite(float(metrics[p+'.'+k])) and float(metrics[p+'.'+k]) >= 0
                   for p in phases for k in phase_components)
        build = builds[role]
        dll = iv['role_trees'][role]['files']['libs/MungeSumstats.so']
        assert dll['sha256'] == build['dll']['sha256']
        assert iv['role_trees'][role]['tree_sha256'] == build['package_tree_sha256']
        assert cfg[role]['commit'] == commits[role] and cfg[role]['dll'] == build['dll']
        assert inp['dll'] == metrics['dll'] == build['dll']['path']
        assert inp['dll_md5'] == metrics['dll_md5'] == dll['md5']
        assert inp['library.0'] == build['library'] and inp['library.1'] == cfg['dependency_libraries'][0]
        base_settings = {'mode':'diagnostic', 'profile':'none', 'allocation_threshold':'65536',
                         'input':case['input']['path'], 'rows':str(n), 'output':str(work),
                         'reference':ref['directory'], 'core_version':f'mss-preservation-{format_version}',
                         'fasta_pin_path':ref['fasta']['path'], 'dll':build['dll']['path'],
                         'dll_md5':dll['md5'], 'library.0':build['library'],
                         'library.1':cfg['dependency_libraries'][0]}
        opts = case['settings']
        base_settings.update({k: str(x) for k,x in opts.items() if k not in {'columns','verify_acquisition_receipt'} and x is not None})
        base_settings.update({f'column.{k}': x for k,x in opts.get('columns',{}).items()})
        if ref['bind_fasta']:
            base_settings['fasta'] = ref['fasta']['path']
        if opts['verify_acquisition_receipt']:
            base_settings.update(acquisition_md5=first['source-artifact']['md5'],
                                 acquisition_bytes=str(first['source-artifact']['bytes']))
        assert inp == base_settings, (case_id, worker, set(inp.items()) ^ set(base_settings.items()))

        raw_ns = ns[case_id, worker]
        assert raw_ns['role'] == role and pin(raw_ns['raw_inventory']) == work/'namespaces.tsv'
        state = table(work/'namespaces.tsv', ['package','version','path'])
        assert len(state) == 71 and state == cert['actual_role_states'][role]
        assert state['MungeSumstats'] == ['1.21.4', str(B/role/'library/MungeSumstats')]
        comparable = {**state,'MungeSumstats':['1.21.4','<selected installation>']}
        assert comparable == cert['common_comparison_state']
        assert comparable == shared.setdefault('state', comparable)
        for name, (version, path) in state.items():
            path = Path(path)
            assert path.resolve(strict=True) == path and path.name == name
            if name != 'MungeSumstats':
                assert name in libraries[str(path.parent)]
            declare = re.search(r'^Version:[ \t]+(\S+)', read(path/'DESCRIPTION').decode('latin-1'), re.M).group(1)
            assert '.'.join(str(int(part)) for part in re.split(r'[.-]',declare)) == version

        files = artifacts[str(work)]
        assert set(files) == {p.name for p in edition.iterdir()} == set(first)
        for name, record in files.items():
            path = edition/name
            assert path.is_file() and not path.is_symlink() and path.stat().st_size == record['bytes']
            if name != 'manifest.rds':
                assert record == first[name]
            if record['bytes'] <= 16384:
                assert digest(read(path)) == record
                if name != 'manifest.rds':
                    assert read(path) == read(D/case_id/order[0]/'edition'/name)
            elif tiny:
                raise AssertionError(('Unexpected large tiny output', path))

        checks = kv(edition/'checksums.tsv', ('filename','md5'))
        assert checks == {name: value['md5'] for name,value in files.items() if name not in {'checksums.tsv','manifest.rds'}}
        assert {k:int(x) for k,x in kv(edition/'audit.tsv', ('status','count')).items()} == v['status_counts']
        if tiny:
            rows = normalized(edition/'normalized.tsv')
            assert rows == gold and len(rows) == n
            assert dict(Counter(row['status'] for row in rows)) == v['status_counts']
            assert read(edition/'source-artifact') == original and read(edition/'source.tsv') == decoded
            tiny_cells += 19*n
        man = gzip.decompress(read(edition/'manifest.rds'))
        assert man.startswith(b'X\n') and man.count(b'created_utc') == 1
        timestamp = rb'(?<=\x00\x00\x00\x17)\d{4}-\d\d-\d\d \d\d:\d\d:\d\d UTC'
        normalized_manifest, count = re.subn(timestamp, b'0'*23, man)
        assert count == 1
        masked.append(normalized_manifest)
        roles[role].append((metrics, process))
    assert len(set(masked)) == 1
    print('SCIENTIFIC', case_id, n, v['status_counts'], 'artifact receipt SHA', sha(D/case_id/'artifact-digests.json'))

    med = {}
    for role, data in roles.items():
        m = {p: {k: statistics.median(float(metrics[p+'.'+k]) for metrics,_ in data)
                 for k in phase_components} for p in phases}
        m.update({summary_key:statistics.median(proc[raw_key] for _,proc in data) for raw_key, summary_key in
                  [('wall_seconds','process_wall_median_seconds'),('user_seconds','process_user_median_seconds'),
                   ('system_seconds','process_system_median_seconds'),('peak_rss_bytes','median_peak_rss_bytes')]})
        m['max_peak_rss_bytes'] = max(proc['peak_rss_bytes'] for _,proc in data)
        m['rss_complete'] = True
        assert m == s['medians'][role]
        med[role] = m
    a,b = med.values()
    assert s['median_peak_rss_nonregression'] == (a['median_peak_rss_bytes'] >= b['median_peak_rss_bytes'])
    pair_operations = {
        'processing':lambda m,p:float(m['processing.elapsed']),
        'wall':lambda m,p:p['wall_seconds'],
        'rss':lambda m,p:p['peak_rss_bytes'],
    }
    pair_results = {key: [f(*a)/f(*b) for a,b in zip(*roles.values())] for key,f in pair_operations.items()}
    assert pair_results['processing'] == s['baseline_over_candidate_by_pair']
    assert a['processing']['elapsed']/b['processing']['elapsed'] == s['ratio_of_median_processing_seconds']
    keys = [('setup', lambda x:x['package_setup']['elapsed']), ('open', lambda x:x['reference_open']['elapsed']),
            ('processing', lambda x:x['processing']['elapsed']), ('wall', lambda x:x['process_wall_median_seconds']),
            ('rss', lambda x:x['median_peak_rss_bytes']), ('maxrss', lambda x:x['max_peak_rss_bytes'])]
    measures = {}
    for label, measure in keys:
        av,bv = measure(a),measure(b)
        measures[label] = [av,bv,av/bv,100*(1-bv/av),av-bv]
    print('MEASURES',case_id,json.dumps(measures))
    print('PAIRS',case_id,json.dumps(pair_results))
    if not tiny:
        print('RAW 1M', {role:[(float(m['processing.elapsed']),p['wall_seconds'],p['peak_rss_bytes'])
                                for m,p in measurements] for role,measurements in roles.items()})
assert tiny_cells == 9348
print('TINY GOLD', tiny_cells, 'field-by-field matches across all tiny workers')

oa = j(O/'admission.json')
oo = observed('oracle', 'cd60999c6f0a3d6c395f6687a4888e742a685bd21d4679a2bbf65ca7e975c6ec', outer=False)
oracle_plan, oracle_run, oracle_terminal = attested_guard(oo, oa['argv'])
assert oa['benchmark_observation'] == {'path':str(T/'mss-key-geometry-benchmark-observation-20261121.json'),
                                       'sha256':'004e670d742c268adf2ca342b90cd41b6e34c85e338de47796dba6b0a93b5775'}
assert oracle_plan['utc'] > bmterm['utc'] and oa['source_commit'] == candidate
assert oa['argv'][1:4] == ['-I','-B','-S'] and sha(oa['argv'][4]) == oa['oracle_sha256']
assert oa['argv'][5:9] == ['--profile',str(pin(oa['profile'])),'--stage','1M']
assert oa['argv'][9:11] == ['--edition',f'mss7={D}/matched-1M/pair-01-candidate/edition']
assert oracle_terminal['minimum_free_bytes'] >= oa['floor_bytes'] and not oracle_terminal['missed_samples']
assert not read(oracle_plan['stderr']) and 'lost across compaction' in oo['outer_tool_terminal']
window = j(pin(oo['input_window']))
oracle = j(pin(oo['oracle']))
assert window['oracle'] == oo['oracle'] and window['equal']
assert window['inputs_before'] == window['inputs_after'] == oa['inputs_before']
assert len(window['inputs_before']) == 5
v = summary['cases']['matched-1M']['verification']
sel = D/'matched-1M/pair-01-candidate'
edition = j(D/'matched-1M/artifact-digests.json')[str(sel)]
assert window['inputs_before'][str(sel/'edition/normalized.tsv')] == edition['normalized.tsv']
assert all(pre['assets'][p] == d for p,d in window['inputs_before'].items() if p != str(sel/'edition/normalized.tsv'))
assert oracle['passed'] and oracle['all_rows_examined'] and oracle['rows_expected'] == 1000000
assert oracle['deltas'] == {} and oracle['unaccounted_deltas'] == 0
assert oracle['source_bytes_sha256'] == edition['source-artifact']['sha256'] == edition['source.tsv']['sha256']
assert oracle['source_copies_exact_each'] == oo['source_copies_exact_each'] == 2
assert oracle['versions']['mss7']['statuses'] == v['status_counts'] == {'match':562542,'swap':437458}
counts = oracle['versions']['mss7']['counts']
assert all(counts[k] == 1000000 for k in ['source_rows','normalized_rows','all_19_fields_correct_rows',
                                          'expected_resolved_rows','observed_resolved_rows'])
assert oo['all_19_fields_correct_rows'] == 1000000 and oracle['reference_inventory']['eligible_panel_records'] == 437458
print('GUARDS', bmo['actual_terminal_tool_chunk'], bmterm['elapsed_seconds'], bmterm['group_samples'],
      bmterm['peak_sampled_group_rss_bytes'], bmterm['peak_disk_allocated_bytes'], bmterm['minimum_free_bytes'],
      oracle_run['pid'], oracle_terminal['elapsed_seconds'], oracle_terminal['group_samples'])
print('PASS independent retained controls')
