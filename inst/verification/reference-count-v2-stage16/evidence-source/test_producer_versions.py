"""Bounded operator provenance tests: owned tiny data and frozen reviewed code."""
import copy
import os
from pathlib import Path
import shutil

import pytest
from test_coordinator import api, fixture, fixture_roots, revise, spec, HERE

PINS = {
    '4': (38909, 'c6327b05066a8a00508c724cd90ac44a0df2ccaa86da43831b43d888dd2b9fb5'),
    '5': (43423, '7a7668c91235ddf8037563ae3141b04e8fd75604c1f24348198a858fb5186147'),
}


def version_fixture(tmp_path, versions=('4', '5')):
    p = fixture(tmp_path)
    p['schema_version'] = 2
    del p['artifacts']['producer']
    for b, version in zip(p['builds'], versions, strict=True):
        path = tmp_path/('producer-v'+version+'.py')
        if not path.exists():
            shutil.copyfile(HERE/'fixture_sources'/('producer_v'+version+'.py'), path)
        b['producer_source'] = spec(path)
        assert (b['producer_source']['size_bytes'], b['producer_source']['sha256']) == PINS[version]
        b['preparation_version'] = 'biodb-reference-stream-'+version
    for i, version in enumerate(versions):
        revise(p, i, 'receipt', lambda v: v.update(preparation_version='biodb-reference-stream-'+version))
    return p


@pytest.mark.parametrize('versions', [('4','4'), ('5','5'), ('4','5')])
def test_schema2_counts_bind_each_exact_source_and_version(tmp_path, versions):
    p = version_fixture(tmp_path, versions)
    before = {str(x): spec(x) for b in p['builds'] for x in Path(b['final_dir']).iterdir()}
    result = api().count_worker(p, tmp_path/'case', expected_publications=fixture_roots(tmp_path))
    assert result['decision'] == 'no_admission'
    assert result['full_reference_admitted'] is False and result['index_builds_executed'] == 0
    for b, r in zip(p['builds'], result['builds'], strict=True):
        assert r['preparation_version'] == b['preparation_version']
        assert {k:r['producer_source'][k] for k in ('path','size_bytes','sha256')} == b['producer_source']
        assert r['counts']['source_rows'] == b['expected_emitted_rows'] == 1
        assert r['producer_source'] in result['controls_before']
        assert r['producer_source'] in result['controls_after']
    assert before == {str(x): spec(x) for b in p['builds'] for x in Path(b['final_dir']).iterdir()}


def refuse_without_counted_access(tmp_path, monkeypatch, p, match=None):
    m = api(); original = Path.open; opened = []
    forbidden = {str(Path(b['final_dir'])/name) for b in p['builds']
                 for name in ('reference.fa','variants.tsv','contig_aliases.tsv')}
    def observe(path, *args, **kwargs):
        opened.append(str(path))
        assert str(path) not in forbidden, 'counted content before refusal'
        return original(path, *args, **kwargs)
    def no_native(*args, **kwargs):
        raise AssertionError('native called before refusal')
    def no_hash(*args, **kwargs):
        raise AssertionError('counted input hash before refusal')
    monkeypatch.setattr(Path, 'open', observe)
    monkeypatch.setattr(m.subprocess, 'run', no_native)
    with pytest.raises(m.Refusal, match=match):
        m.count_worker(p, tmp_path/'case', input_pinner=no_hash, expected_publications=fixture_roots(tmp_path))
    assert not (tmp_path/'case'/'worker.json').exists()
    return opened


@pytest.mark.parametrize('location', ['receipt','expected'])
@pytest.mark.parametrize('bad', [None, 'biodb-reference-stream-0', 5, True, 'opposite', 'missing'])
def test_wrong_missing_or_unknown_version_refuses_before_inputs(tmp_path, monkeypatch, location, bad):
    p = version_fixture(tmp_path)
    bad = 'biodb-reference-stream-5' if bad == 'opposite' else bad
    if location == 'receipt':
        revise(p, 0, 'receipt', lambda v: v.pop('preparation_version') if bad == 'missing' else v.update(preparation_version=bad))
    elif bad == 'missing':
        p['builds'][0].pop('preparation_version')
    else:
        p['builds'][0]['preparation_version'] = bad
    refuse_without_counted_access(tmp_path, monkeypatch, p)


@pytest.mark.parametrize('index', [0,1])
def test_receipt_version_swapped_refuses(tmp_path, monkeypatch, index):
    p = version_fixture(tmp_path)
    revise(p, index, 'receipt', lambda v: v.update(preparation_version='biodb-reference-stream-'+('5' if index==0 else '4')))
    refuse_without_counted_access(tmp_path, monkeypatch, p)


@pytest.mark.parametrize('failure', ['missing37','missing38','swapped','sha','size','bool','type',
                                    'extra','global_only','redundant','missing_common'])
def test_bad_source_authority_refuses_before_inputs(tmp_path, monkeypatch, failure):
    p = version_fixture(tmp_path)
    b = p['builds'][0]
    if failure.startswith('missing') and failure != 'missing_common':
        p['builds'][int(failure=='missing38')].pop('producer_source')
    elif failure == 'swapped':
        b['producer_source'], p['builds'][1]['producer_source'] = p['builds'][1]['producer_source'], b['producer_source']
    elif failure in ('sha','size','bool'):
        b['producer_source'][{'sha':'sha256','size':'size_bytes','bool':'size_bytes'}[failure]] = {'sha':'a'*64,'size':43423,'bool':True}[failure]
    elif failure == 'type': b['producer_source'] = 'wrong'
    elif failure == 'extra': b['producer_source']['version'] = '4'
    elif failure in ('global_only','redundant'):
        p['artifacts']['producer'] = copy.deepcopy(b['producer_source'])
        if failure == 'global_only':
            for item in p['builds']: item.pop('producer_source')
    elif failure == 'missing_common': p['artifacts'].pop('header')
    refuse_without_counted_access(tmp_path, monkeypatch, p)


@pytest.mark.parametrize('failure', ['same_size','wrong_size','huge'])
def test_actual_source_mismatch_refuses_without_counted_access(tmp_path, monkeypatch, failure):
    p = version_fixture(tmp_path); path = Path(p['builds'][0]['producer_source']['path'])
    if failure == 'same_size': path.write_bytes(b'x'*38909)
    else:
        with path.open('wb') as handle: handle.truncate(38910 if failure=='wrong_size' else 32*1024**2)
    original = Path.open
    def bounded(target,*args,**kwargs):
        assert failure=='same_size' or target!=path, 'wrong observed source size opened'
        return original(target,*args,**kwargs)
    monkeypatch.setattr(Path, 'open', bounded)
    refuse_without_counted_access(tmp_path, monkeypatch, p)


@pytest.mark.parametrize('failure', ['inside','alias_data','hardlink','symlink'])
def test_source_metadata_placement_refuses_before_content(tmp_path, monkeypatch, failure):
    p = version_fixture(tmp_path); b = p['builds'][0]
    old = Path(b['producer_source']['path']); target = tmp_path/'bad-source.py'
    if failure=='inside': target=Path(b['final_dir'])/'source.py';shutil.copyfile(old,target)
    elif failure=='alias_data': os.link(Path(b['final_dir'])/'variants.tsv', target)
    elif failure=='hardlink': os.link(old,target)
    else: target.symlink_to(old)
    b['producer_source']['path']=str(target)
    opened=[]
    def no_open(path,*args,**kwargs):
        opened.append(str(path));raise AssertionError('metadata refusal opened content')
    monkeypatch.setattr(Path,'open',no_open)
    m=api()
    with pytest.raises(m.Refusal):m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
    assert opened==[]


@pytest.mark.parametrize('failure', ['running','failed','wrong_command','exit','oom'])
def test_both_guards_precede_new_source_and_all_seal_content(tmp_path, monkeypatch, failure):
    p=version_fixture(tmp_path)
    changes={'running':{'status':'running'},'failed':{'status':'failed'},'wrong_command':{'command':['wrong']},'exit':{'exit_code':1},'oom':{'natural_oom':True}}
    revise(p,1,'producer_guard',lambda v:v.update(changes[failure]))
    guards={b['producer_guard']['path'] for b in p['builds']}; original=Path.open;opened=[]
    def observe(path,*args,**kwargs):
        opened.append(str(path));assert str(path) in guards,'content before both guards'
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    m=api()
    with pytest.raises(m.Refusal):m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
    assert set(opened)==guards


def test_source_mutation_during_worker_prevents_success(tmp_path, monkeypatch):
    p=version_fixture(tmp_path);m=api();original=m.subprocess.run; changed=False
    def mutate(*args,**kwargs):
        nonlocal changed
        result=original(*args,**kwargs)
        if not changed:
            changed=True;path=Path(p['builds'][0]['producer_source']['path'])
            path.write_bytes(b'x'*38909)
        return result
    monkeypatch.setattr(m.subprocess,'run',mutate)
    with pytest.raises(m.Refusal):m.count_worker(p,tmp_path/'case',expected_publications=fixture_roots(tmp_path))
    assert changed and not (tmp_path/'case'/'worker.json').exists()


def test_v1_still_rejects_v5_before_inputs(tmp_path, monkeypatch):
    p=fixture(tmp_path)
    p['builds'][0]['preparation_version']='biodb-reference-stream-5'
    revise(p,0,'receipt',lambda v:v.update(preparation_version='biodb-reference-stream-5'))
    refuse_without_counted_access(tmp_path,monkeypatch,p)


@pytest.mark.parametrize('target', ['source_plan','producer_source'])
def test_schema2_repeated_controls_check_bounds_before_open(tmp_path,monkeypatch,target):
    p=version_fixture(tmp_path);m=api();original=m.validate_seals;original_open=Path.open;changed=False
    path=Path(p['builds'][0][target]['path'])
    def mutate(*args,**kwargs):
        nonlocal changed
        result=original(*args,**kwargs)
        with original_open(path,'wb') as handle:handle.truncate(16*1024**2+1)
        changed=True
        return result
    def observe(candidate,*args,**kwargs):
        assert not (changed and candidate==path),'oversize repeated control opened'
        return original_open(candidate,*args,**kwargs)
    monkeypatch.setattr(m,'validate_seals',mutate);monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.count_worker(p,tmp_path/'case',expected_publications=fixture_roots(tmp_path))
    assert changed and not (tmp_path/'case'/'worker.json').exists()


def test_both_complete_seals_precede_source_content_and_counted_inputs(tmp_path,monkeypatch):
    p=version_fixture(tmp_path)
    revise(p,1,'receipt',lambda v:v.update(alias_history_status='not_supplied'))
    forbidden={b['producer_source']['path'] for b in p['builds']}
    original=Path.open
    def observe(path,*args,**kwargs):
        assert str(path) not in forbidden,'producer code opened before both complete seals'
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    refuse_without_counted_access(tmp_path,monkeypatch,p)


def test_trusted_command_can_import_module_without_source_path_in_argv(tmp_path):
    p=version_fixture(tmp_path)
    # Fixture's trusted command contains plan/publication, no module path argument.
    assert all(b['producer_source']['path'] not in b['producer_command'] for b in p['builds'])
    result=api().count_worker(p,tmp_path/'case',expected_publications=fixture_roots(tmp_path))
    assert result['decision']=='no_admission'


@pytest.mark.parametrize('failure',['expected_large','expected_bool','observed_large'])
def test_schema2_json_expected_observed_bounds_before_content(tmp_path,monkeypatch,failure):
    p=version_fixture(tmp_path);b=p['builds'][0];source=Path(b['producer_guard']['path'])
    if failure=='observed_large':
        with source.open('wb') as handle:handle.truncate(16*1024**2+1)
    else:b['producer_guard']['size_bytes']=True if failure=='expected_bool' else 16*1024**2+1
    def no_open(*args,**kwargs):raise AssertionError('overbound JSON opened')
    monkeypatch.setattr(Path,'open',no_open)
    m=api()
    with pytest.raises(m.Refusal):m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
