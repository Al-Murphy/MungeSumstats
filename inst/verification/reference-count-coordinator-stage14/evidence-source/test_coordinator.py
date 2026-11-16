"""Tiny owned fixtures only: no production paths or preparation output reads."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
MODULE = HERE / 'sealed_count.py'


def api():
    assert MODULE.exists(), 'sealed coordinator implementation missing'
    spec = importlib.util.spec_from_file_location('sealed_count', MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def spec(path):
    return {'path': str(path), 'size_bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def write(path, value):
    path.write_text(json.dumps(value) + '\n')
    return spec(path)


def fixture(tmp_path):
    """Synthetic authenticated local context; never promoted to full scope."""
    roles = ['coordinator', 'counter_source', 'native', 'header', 'guard', 'harness', 'producer']
    pins = {}
    actual = {'coordinator': MODULE, 'counter_source': HERE / 'count_reference.cpp',
              'native': HERE / 'count_reference'}
    for role in roles:
        p = actual.get(role, tmp_path / (role + '.txt'))
        if role not in actual:
            p.write_text(role)
        pins[role] = spec(p)
    builds = []
    for build in ['GRCh37', 'GRCh38']:
        final = tmp_path / build
        final.mkdir()
        (final / 'reference.fa').write_bytes(b'>1\r\nAAAAAAAA\r\n')
        (final / 'variants.tsv').write_bytes(b'chr\tpos\tref\talt\tid\taliases\n1\t1\tA\tC,G,AA\trs1\trs1;old;old\n')
        (final / 'contig_aliases.tsv').write_text('alias\tcanonical_contig\nchr1\t1\n')
        sources = {role: {'sha256': 'a'*64, 'size_bytes': 1, 'build': build, 'release': '157',
                          'completeness': 'complete_source'}
                   for role in ['fasta', 'vcf', 'assembly_report', 'merged', 'withdrawn']}
        plan = {'build': build, 'release': 'tiny-v1', 'scope': 'tiny_fixture',
                'assembly': 'GCF-test', 'vcf_headers': ['##reference='+build], 'sources': sources}
        plan_pin = write(tmp_path / (build + '-plan.json'), plan)
        outputs = {name: {k: v for k, v in spec(final/name).items() if k != 'path'}
                   for name in ['reference.fa', 'variants.tsv', 'contig_aliases.tsv']}
        receipt = dict(plan, preparation_version='biodb-reference-stream-4',
                       reference_format='mss-preservation-2', alias_history_status='complete_supplied_history',
                       reference_fasta=str(final/'reference.fa'), contig_aliases=str(final/'contig_aliases.tsv'),
                       counts={'emitted_rows': 1}, outputs=outputs)
        command = ['tiny-producer', '--manifest', plan_pin['path'], '--output-dir', str(final)]
        guard = {'status': 'success', 'exit_code': 0, 'reason': None, 'natural_oom': False,
                 'command': command}
        builds.append({'build': build, 'release': 'tiny-v1', 'final_dir': str(final),
                       'source_plan': plan_pin, 'receipt': write(final/'receipts.json', receipt),
                       'producer_guard': write(tmp_path/(build+'-guard.json'), guard),
                       'producer_command': command, 'preparation_version': 'biodb-reference-stream-4',
                       'expected_emitted_rows': 1})
    return {'schema_version': 1, 'context': 'tiny_fixture', 'format_version': 3,
            'artifacts': pins, 'builds': builds, 'manifest_upper_bytes': 1048576,
            'manifest_bound_basis': 'Explicit operator planning allowance; no RDS build qualification',
            'allocation_allowance_bytes': 4096, 'concurrent_reservation_bytes': 0,
            'limits': {'tree_memory_bytes': 8*1024**3, 'minimum_free_bytes': 160*1024**3,
                       'owned_case_bytes': 16*1024**2, 'deadline_seconds': 60}}


def fixture_roots(tmp_path, published=False):
    """Independent paths fixed by fixture construction, outside profile JSON."""
    suffix = '-published' if published else ''
    return {build: str(tmp_path/(build+suffix)) for build in ('GRCh37', 'GRCh38')}


def revise(profile, build_index, kind, mutate):
    item = profile['builds'][build_index]
    path = Path(item[kind]['path'])
    value = json.loads(path.read_text())
    mutate(value)
    item[kind] = write(path, value)  # Simulates reviewed context pinning BAD semantic metadata.


def test_both_seals_before_input_open(tmp_path):
    m = api(); p = fixture(tmp_path)
    revise(p, 1, 'producer_guard', lambda v: v.update(status='running'))
    opened = []
    with pytest.raises(m.Refusal):
        m.count_worker(p, tmp_path/'case', input_pinner=lambda x: opened.append(x), expected_publications=fixture_roots(tmp_path))
    assert opened == []


@pytest.mark.parametrize('kind,change', [
    ('producer_guard', lambda v: v.update(exit_code=1)),
    ('producer_guard', lambda v: v.update(command=['wrong'])),
    ('producer_guard', lambda v: v.update(natural_oom=True)),
    ('receipt', lambda v: v.update(build='GRCh38')),
    ('receipt', lambda v: v.update(release='wrong')),
    ('receipt', lambda v: v.update(preparation_version='v0')),
    ('receipt', lambda v: v.update(alias_history_status='not_supplied')),
    ('receipt', lambda v: v['sources']['merged'].update(completeness='prefix')),
    ('receipt', lambda v: v['outputs'].pop('variants.tsv')),
    ('receipt', lambda v: v['counts'].update(emitted_rows=True)),
    ('source_plan', lambda v: v.update(release='altered')),
])
def test_semantic_seal_refusals(tmp_path, kind, change):
    m = api(); p = fixture(tmp_path)
    revise(p, 0, kind, change)
    with pytest.raises(m.Refusal):
        m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def test_unreviewed_digest_change_refused(tmp_path):
    m = api(); p = fixture(tmp_path)
    Path(p['builds'][1]['receipt']['path']).write_text('{}')
    with pytest.raises(m.Refusal): m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def test_symlink_final_refused(tmp_path):
    m = api(); p = fixture(tmp_path)
    link = tmp_path/'linked'; link.symlink_to(p['builds'][0]['final_dir'])
    p['builds'][0]['final_dir'] = str(link)
    with pytest.raises(m.Refusal): m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def test_full_context_requires_separate_admission(tmp_path):
    m = api(); p = fixture(tmp_path); p['context'] = 'full_reference'
    with pytest.raises(m.Refusal): m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def counts(**kw):
    result = dict(source_bytes=100, records_bytes=70, raw_sort_bytes=120, keys_bytes=90,
                  final_key_rows=10, rejected_bytes=20, fai_upper_bytes=86)
    result.update(kw)
    return result


def test_exact_arithmetic_and_both_orders():
    m = api()
    a = m.logical_size(counts(), 10, 100)
    assert a == {'small_files_upper_bytes': 215, 'retained_bytes': 575, 'phase_bytes': 595}
    b = m.logical_size(counts(raw_sort_bytes=600), 0, 0)
    assert b == {'small_files_upper_bytes': 105, 'retained_bytes': 465, 'phase_bytes': 1395}
    plan = m.capacity([a, b], 170*1024**3 + 1980, 10, 0)
    assert plan['order_peak_bytes'] == {'GRCh37_then_GRCh38': 1970, 'GRCh38_then_GRCh37': 1395}
    assert plan['orders_fit_with_reserves'] == {'GRCh37_then_GRCh38': True, 'GRCh38_then_GRCh37': True}
    assert not m.capacity([a, b], 170*1024**3+1979, 10, 0)['orders_fit_with_reserves']['GRCh37_then_GRCh38']
    assert plan['full_reference_admitted'] is False and plan['index_builds_executed'] == 0


@pytest.mark.parametrize('value', [-1, True, 1.0, '2', 2**64])
def test_exact_integer_refusal(value):
    m = api()
    with pytest.raises(m.Refusal): m.logical_size(counts(final_key_rows=value), 0, 0)


def test_checked_sum_overflow():
    m = api()
    with pytest.raises(m.Refusal): m.logical_size(counts(source_bytes=2**63-1, records_bytes=2**63-1), 0, 0)


def test_arithmetic_random_oracle():
    m = api(); rng = random.Random(17)
    for _ in range(200):
        t,r,u,k,n,j,f,c,d = [rng.randrange(10000) for _ in range(9)]
        actual = m.logical_size(counts(source_bytes=t, records_bytes=r, raw_sort_bytes=u,
                                      keys_bytes=k, final_key_rows=n, rejected_bytes=j, fai_upper_bytes=f), c, d)
        small=f+c+d+19
        assert actual['retained_bytes'] == t+r+k+8*n+j+small
        assert actual['phase_bytes'] == max(t+r+k+8*n+j+small,r+2*u+j+small,r+u+k+8*n+j+small)


def test_alias_bound_quoted_unicode(tmp_path):
    m = api(); p = tmp_path/'aliases.tsv'
    p.write_text('alias\tcanonical_contig\n"a""b"\t染色体\n')
    assert m.alias_upper(p) == len('"alias"\t"canonical_contig"\r\n"a""b"\t"染色体"\r\n'.encode())


def test_actual_count_worker_readonly_exclusive(tmp_path):
    m = api(); p = fixture(tmp_path); case = tmp_path/'case'
    before = {str(f): spec(f) for b in p['builds'] for f in Path(b['final_dir']).iterdir()}
    result = m.count_worker(p, case, expected_publications=fixture_roots(tmp_path))
    assert [x['counts']['source_rows'] for x in result['builds']] == [1, 1]
    assert result['full_reference_admitted'] is False and result['index_builds_executed'] == 0
    assert before == {str(f): spec(f) for b in p['builds'] for f in Path(b['final_dir']).iterdir()}
    assert sorted(x.name for x in case.iterdir()) == ['GRCh37.stderr', 'GRCh37.stdout', 'GRCh38.stderr', 'GRCh38.stdout', 'worker.json']
    with pytest.raises(FileExistsError): m.count_worker(p, case, expected_publications=fixture_roots(tmp_path))


@pytest.mark.parametrize('failure', ['native_error', 'stderr', 'missing', 'wrong_rows', 'invalid_number', 'mutation'])
def test_native_failures_and_input_mutation(tmp_path, failure):
    m = api(); p = fixture(tmp_path)
    script = tmp_path/'fake-native'
    real = str(HERE/'count_reference')
    body = '#!'+sys.executable+'\nimport json,subprocess,sys\n'
    if failure == 'native_error': body += 'sys.exit(2)\n'
    else:
        body += 'v=json.loads(subprocess.check_output(['+repr(real)+',*sys.argv[1:]]))\n'
        body += {'stderr':'print("warning",file=sys.stderr)\n', 'missing':'v.pop("keys_bytes")\n',
                 'wrong_rows':'v["source_rows"]=2\n','invalid_number':'v["keys_bytes"]=True\n',
                 'mutation':'with open(sys.argv[2],"ab") as f: f.write(b"\\n")\n'}[failure]
        body += 'print(json.dumps(v))\n'
    script.write_text(body); script.chmod(0o700); p['artifacts']['native'] = spec(script)
    with pytest.raises(m.Refusal): m.count_worker(p, tmp_path/'case', expected_publications=fixture_roots(tmp_path))
    assert not (tmp_path/'case'/'worker.json').exists()


def test_guard_censor_never_success(tmp_path):
    m = api()
    for status in ['censored', 'failed']:
        with pytest.raises(m.Refusal):
            m.check_guard({'status':status,'exit_code':0,'reason':None,'natural_oom':False})


def test_entire_source_plan_policy_bound(tmp_path):
    m = api(); p = fixture(tmp_path)
    revise(p, 0, 'source_plan', lambda v: v.update(ref_mismatch_policy='quarantine_compatible_ambiguity'))
    with pytest.raises(m.Refusal): m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def test_pending_path_refused(tmp_path):
    m = api(); p = fixture(tmp_path)
    final = Path(p['builds'][0]['final_dir']); pending = final.with_name('.GRCh37-pending-test')
    final.rename(pending); p['builds'][0]['final_dir'] = str(pending)
    with pytest.raises(m.Refusal): m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def test_native_addressability_boundary():
    m = api()
    with pytest.raises(m.Refusal): m.logical_size(counts(final_key_rows=m.ADDRESS_MAX//8+1), 0, 0)


def test_deadline_bool_and_relaxed_guard_refused(tmp_path):
    m = api(); p = fixture(tmp_path)
    for key, value in [('deadline_seconds', True), ('minimum_free_bytes', 1), ('tree_memory_bytes', 2**64)]:
        bad = copy.deepcopy(p); bad['limits'][key] = value
        with pytest.raises(m.Refusal): m.validate_seals(bad, expected_publications=fixture_roots(tmp_path))


def test_same_size_output_digest_change_refused(tmp_path):
    m = api(); p = fixture(tmp_path)
    path = Path(p['builds'][0]['final_dir'])/'variants.tsv'
    path.write_bytes(path.read_bytes().replace(b'old', b'new'))
    with pytest.raises(m.Refusal): m.count_worker(p, tmp_path/'case', expected_publications=fixture_roots(tmp_path))


def test_pinned_json_rejects_duplicate_keys(tmp_path):
    m = api(); path = tmp_path/'duplicate.json'; path.write_text('{"a":1,"a":2}')
    with pytest.raises(m.Refusal): m.pinned_json(spec(path))


def test_cli_cannot_bypass_monitor(tmp_path):
    m = api(); p = fixture(tmp_path)
    profile = tmp_path/'profile.json'; pp = write(profile, p)
    result = subprocess.run([sys.executable, str(MODULE), '--worker', '--profile', str(profile),
                             '--expected-profile-sha256', pp['sha256'], '--case', str(tmp_path/'unmonitored')], capture_output=True)
    assert result.returncode != 0
    assert not (tmp_path/'unmonitored').exists()


def test_full_admission_requires_same_frozen_source(tmp_path):
    m = api(); p = fixture(tmp_path)
    p.update(context='full_reference', coordinator_source_commit='a'*40)
    admission = dict(action='sealed-full-count-only', independent_review='GREEN',
                     source_commit='b'*40, profile_sha256='c'*64)
    opened=[]
    with pytest.raises(m.Refusal, match='source'):
        m.count_worker(p,tmp_path/'case',admission,input_pinner=lambda x:opened.append(x), expected_publications=fixture_roots(tmp_path))
    assert opened==[]


def test_invalid_native_json_refusal(tmp_path):
    m=api();p=fixture(tmp_path);script=tmp_path/'invalid-json'
    script.write_text('#!'+sys.executable+'\nprint("{truncated")\n');script.chmod(0o700)
    p['artifacts']['native']=spec(script)
    with pytest.raises(m.Refusal):m.count_worker(p,tmp_path/'case', expected_publications=fixture_roots(tmp_path))


def test_misplaced_artifact_never_opens_input_before_both_seals(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path)
    target=Path(p['builds'][0]['final_dir'])/'variants.tsv'
    p['artifacts']['native']=spec(target)
    revise(p,1,'producer_guard',lambda v:v.update(status='running'))
    original=Path.open
    def checked(path,*args,**kwargs):
        assert path!=target,'prepared input opened before both seals'
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',checked)
    with pytest.raises(m.Refusal):m.validate_seals(p, expected_publications=fixture_roots(tmp_path))


def test_hash_refuses_mutation_during_read(tmp_path,monkeypatch):
    m=api();target=tmp_path/'changing';target.write_bytes(b'a'*100)
    original=Path.open
    class ChangingReader:
        def __init__(self,handle):self.handle=handle;self.changed=False
        def __enter__(self):return self
        def __exit__(self,*args):self.handle.close()
        def fileno(self):return self.handle.fileno()
        def read(self,n):
            value=self.handle.read(n)
            if not self.changed:
                self.changed=True
                with original(target,'ab') as h:h.write(b'b')
            return value
    def checked(path,*args,**kwargs):
        handle=original(path,*args,**kwargs)
        return ChangingReader(handle) if path==target and args==('rb',) else handle
    monkeypatch.setattr(Path,'open',checked)
    with pytest.raises(m.Refusal,match='changed'):m.pin(target)


def test_final_receipt_allowance_enforced_before_writing(tmp_path):
    m=api();target=tmp_path/'summary.json'
    with pytest.raises(m.Refusal):m.write_once(target,{'value':'x'*100},limit_bytes=10)
    assert not target.exists()
