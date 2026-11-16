"""Focused regressions for the four independently demonstrated coordinator gaps."""
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from test_coordinator import api, fixture, spec, write, HERE

CAP=16*1024**2


def roots(path, published=False):
    suffix='-published' if published else ''
    return {build:str(path/(build+suffix)) for build in ('GRCh37','GRCh38')}


def tree_names(path):
    return sorted(str(p.relative_to(path)) for p in path.rglob('*'))


@pytest.mark.parametrize('kind',['case37','case38','profile37','admission38','case_symlink',
                                  'profile_symlink','admission_hardlink','missing_context'])
def test_external_placement_refuses_before_any_control_open_or_mkdir(tmp_path,monkeypatch,kind):
    m=api();p=fixture(tmp_path);expected=roots(tmp_path)
    profile=tmp_path/'profile.json';admission=None;case=tmp_path/'case'
    if kind.startswith('case') and kind!='case_symlink':case=tmp_path/('GRCh37' if kind=='case37' else 'GRCh38')/'case'
    if kind=='profile37':profile=tmp_path/'GRCh37'/'profile.json'
    if kind=='admission38':admission=tmp_path/'GRCh38'/'admission.json';write(admission,{})
    if kind in ('case_symlink','profile_symlink'):
        alias=tmp_path/'source-alias';alias.symlink_to(tmp_path/'GRCh37')
        if kind=='case_symlink':case=alias/'case'
        else:profile=alias/'profile.json'
    if kind=='admission_hardlink':
        inside=tmp_path/'GRCh38'/'admission.json';write(inside,{})
        admission=tmp_path/'admission-link.json';os.link(inside,admission)
    pp=write(profile,p)
    before={b:tree_names(Path(d)) for b,d in expected.items()}
    opened=[];original=Path.open
    def observe(path,*args,**kwargs):
        opened.append(str(path));raise AssertionError('control opened before placement refusal: '+str(path))
    monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):
        m.launch(profile,pp['sha256'],case,admission,'a'*64 if admission else None,
                 expected_publications=None if kind=='missing_context' else expected)
    assert opened==[] and not case.exists()
    assert before=={b:tree_names(Path(d)) for b,d in expected.items()}


@pytest.mark.parametrize('entry',['worker','seals'])
def test_internal_entry_requires_independent_context_before_open(tmp_path,monkeypatch,entry):
    m=api();p=fixture(tmp_path)
    def forbidden(*args,**kwargs):raise AssertionError('missing context opened content')
    monkeypatch.setattr(Path,'open',forbidden)
    with pytest.raises(m.Refusal):
        if entry=='worker':m.count_worker(p,tmp_path/'case')
        else:m.validate_seals(p)
    assert not (tmp_path/'case').exists()


def test_internal_worker_case_placement_before_controls(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path);expected=roots(tmp_path);before=tree_names(tmp_path/'GRCh37')
    def forbidden(*args,**kwargs):raise AssertionError('misplaced worker opened a control')
    monkeypatch.setattr(Path,'open',forbidden)
    with pytest.raises(m.Refusal):m.count_worker(p,tmp_path/'GRCh37'/'case',expected_publications=expected)
    assert tree_names(tmp_path/'GRCh37')==before


def test_parsed_roots_must_equal_independent_context_before_seal_open(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path);expected=roots(tmp_path);expected['GRCh37'],expected['GRCh38']=expected['GRCh38'],expected['GRCh37']
    profile=tmp_path/'profile.json';pp=write(profile,p);original=Path.open;opened=[]
    def observe(path,*args,**kwargs):
        opened.append(str(path));assert path==profile,'seal/input opened before expected-root equality'
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.load_profile(profile,pp['sha256'],expected_publications=expected)
    assert opened and set(opened)=={str(profile)}


@pytest.mark.parametrize('kind',['expected_large','observed_large','expected_bool'])
def test_control_expected_and_observed_bounds_precede_any_content_open(tmp_path,monkeypatch,kind):
    m=api();path=tmp_path/'bounded.json'
    if kind=='observed_large':
        with path.open('wb') as handle:handle.truncate(CAP+1)
        pin={'path':str(path),'size_bytes':2,'sha256':'a'*64}
    else:
        path.write_text('{}');pin=spec(path)
        pin['size_bytes']=True if kind=='expected_bool' else CAP+1
    opens=[]
    def forbidden(path,*args,**kwargs):opens.append(str(path));raise AssertionError('oversize control content opened')
    monkeypatch.setattr(Path,'open',forbidden)
    with pytest.raises(m.Refusal):m.pinned_json(pin)
    assert opens==[]


def test_exact_control_cap_is_positive(tmp_path):
    m=api();path=tmp_path/'exact-cap.json'
    path.write_bytes(b'"'+b'a'*(CAP-2)+b'"')
    result=m.pinned_json(spec(path))
    assert len(result)==CAP-2


def test_root_context_cli_is_required(tmp_path):
    result=subprocess.run([sys.executable,str(HERE/'sealed_count.py'),'--profile',str(tmp_path/'not-opened'),
                           '--expected-profile-sha256','a'*64,'--case',str(tmp_path/'case')],capture_output=True)
    assert result.returncode!=0 and not (tmp_path/'case').exists()
    assert b'expected-grch37-publication' in result.stderr


def native_positive(tmp_path,rows=None):
    fasta=tmp_path/'native.fa';fasta.write_bytes(b'>1\nAAAAAAAA\n')
    source=tmp_path/'native.tsv'
    rows=['1\t1\tA\tC,G,AA\trs1\trs1;old;old'] if rows is None else rows
    source.write_text('chr\tpos\tref\talt\tid\taliases\n'+'\n'.join(rows)+ ('\n' if rows else ''))
    run=subprocess.run([str(HERE/'count_reference'),str(fasta),str(source),'3'],capture_output=True)
    assert run.returncode==0 and not run.stderr
    value=json.loads(run.stdout)
    seal={'expected_emitted_rows':len(rows),'inputs':[{},spec(source)]}
    return value,seal


BAD_COUNTS=[
    {'records_bytes':0}, {'raw_alt_tokens':0}, {'raw_key_rows':0,'raw_sort_bytes':0,'final_key_rows':0,'keys_bytes':0,'offsets_bytes':0},
    {'rejected_bytes':0}, {'fai_upper_bytes':0}, {'max_record_bytes':0}, {'max_record_bytes':1000000},
    {'max_key_row_bytes':0}, {'max_key_row_bytes':1000000}, {'max_source_line_bytes':0},
    {'max_source_line_bytes':1000000}, {'max_raw_keys_per_record':0}, {'max_raw_keys_per_record':1000000},
    {'raw_alt_tokens':4}, {'rejected_alts':1}, {'rejected_bytes':12},
    {'records_bytes':1}, {'keys_bytes':1}, {'raw_sort_bytes':1},
]


@pytest.mark.parametrize('change',BAD_COUNTS)
def test_malformed_native_neighbors_refused_before_planning(tmp_path,change):
    m=api();value,seal=native_positive(tmp_path);value.update(change)
    with pytest.raises(m.Refusal):m.validate_counts(value,seal)


@pytest.mark.parametrize('rows',[
    [], ['1\tbad\tA\tC\trs1\t.'], ['1\tbad\tA\t\t.\t.'],
    ['1\t1\tA\t\t.\t.'], ['1\t1\tA\tN,<DEL>,*,A,\trs1\t.'],
    ['1\t1\tA\tC,G,AA\trs1\trs1;old;old'],
    ['1\tbad\tA\tC,G\trs1\t.','1\t1\tA\tC\t.\t.'],
    ['1\t1\tA\tC\t.\t.'], ['1\t1\tA\tC\trs1\t.'],
    ['\t1\t\t\t\t'], ['1\t1\tA\tC\t\t'],
    ['1\t1\tA\tC\t識別子\t旧名;旧名'],
])
def test_valid_empty_rejected_alt_and_alias_boundaries_retained(tmp_path,rows):
    m=api();value,seal=native_positive(tmp_path,rows)
    m.validate_counts(value,seal)
    assert value['source_rows']==len(rows)


def test_empty_maximum_and_rejection_boundaries(tmp_path):
    m=api();value,seal=native_positive(tmp_path,[])
    for key in ['record_rows','raw_key_rows','final_key_rows','max_record_bytes','max_key_row_bytes',
                'max_raw_keys_per_record','max_source_line_bytes','raw_alt_tokens','rejected_alts']:
        bad=dict(value);bad[key]=1
        with pytest.raises(m.Refusal):m.validate_counts(bad,seal)


def test_impossible_injected_response_leaves_no_worker_success(tmp_path):
    m=api();p=fixture(tmp_path);native=tmp_path/'corrupt-native'
    native.write_text('#!'+sys.executable+'\nimport json,subprocess,sys\n'
                      'v=json.loads(subprocess.check_output(['+repr(str(HERE/'count_reference'))+',*sys.argv[1:]]))\n'
                      'for k in ["records_bytes","keys_bytes","raw_sort_bytes","raw_key_rows","final_key_rows","offsets_bytes","rejected_bytes","fai_upper_bytes","raw_alt_tokens"]:v[k]=0\n'
                      'print(json.dumps(v))\n')
    native.chmod(0o700);p['artifacts']['native']=spec(native)
    with pytest.raises(m.Refusal):m.count_worker(p,tmp_path/'case',expected_publications=roots(tmp_path))
    assert not (tmp_path/'case'/'worker.json').exists()


def test_repeated_control_hash_rechecks_observed_bound(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path);expected=roots(tmp_path)
    control=Path(p['builds'][0]['source_plan']['path']);original_open=Path.open
    original_validate=m.validate_seals;changed=False;opens=[]
    def mutate_after_seals(*args,**kwargs):
        nonlocal changed
        result=original_validate(*args,**kwargs)
        with original_open(control,'wb') as handle:handle.truncate(CAP+1)
        changed=True
        return result
    def observe(path,*args,**kwargs):
        if changed and path==control:
            opens.append(str(path));raise AssertionError('repeat hash opened an oversized control')
        return original_open(path,*args,**kwargs)
    monkeypatch.setattr(m,'validate_seals',mutate_after_seals);monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.count_worker(p,tmp_path/'case',expected_publications=expected)
    assert opens==[] and not (tmp_path/'case'/'worker.json').exists()


def test_external_rehash_rechecks_observed_bound_before_hash(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path);expected=roots(tmp_path)
    path=tmp_path/'profile.json';pp=write(path,p);original_open=Path.open
    original_load=m.load_profile;changed=False;opens=[]
    def mutate_after_load(*args,**kwargs):
        nonlocal changed
        result=original_load(*args,**kwargs)
        with original_open(path,'wb') as handle:handle.truncate(CAP+1)
        changed=True
        return result
    def observe(target,*args,**kwargs):
        if changed and target==path:
            opens.append(str(target));raise AssertionError('external rehash opened an oversized control')
        return original_open(target,*args,**kwargs)
    monkeypatch.setattr(m,'load_profile',mutate_after_load);monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.launch(path,pp['sha256'],tmp_path/'case',expected_publications=expected)
    assert opens==[] and not (tmp_path/'case'/'summary.json').exists()


def test_minimal_empty_header_and_crlf_native_positives(tmp_path):
    m=api();value,seal=native_positive(tmp_path,[])
    source=Path(seal['inputs'][1]['path'])
    for data in [b'chr\tpos\tref\talt\tid\taliases',
                 b'chr\tpos\tref\talt\tid\taliases\r\n1\t1\tA\tC\t\t\r\n']:
        source.write_bytes(data)
        result=subprocess.run([str(HERE/'count_reference'),str(tmp_path/'native.fa'),str(source),'3'],capture_output=True)
        assert result.returncode==0 and not result.stderr
        actual=json.loads(result.stdout)
        seal={'expected_emitted_rows':actual['source_rows'],'inputs':[{},spec(source)]}
        m.validate_counts(actual,seal)


def test_expected_root_symlink_alias_refused_before_control_open(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path);expected=roots(tmp_path)
    alias=tmp_path/'root-alias';alias.symlink_to(tmp_path/'GRCh37');expected['GRCh37']=str(alias)
    path=tmp_path/'profile.json';pp=write(path,p)
    def forbidden(*args,**kwargs):raise AssertionError('aliased expected context opened profile')
    monkeypatch.setattr(Path,'open',forbidden)
    with pytest.raises(m.Refusal):m.launch(path,pp['sha256'],tmp_path/'case',expected_publications=expected)
    assert not (tmp_path/'case').exists()


@pytest.mark.parametrize('maximum,total',[('max_record_bytes','records_bytes'),('max_key_row_bytes','keys_bytes')])
def test_nonempty_neighbor_records_require_their_own_minimum_bytes(tmp_path,maximum,total):
    m=api();value,seal=native_positive(tmp_path);value[maximum]=value[total]
    with pytest.raises(m.Refusal):m.validate_counts(value,seal)


def test_other_source_rows_require_their_own_framing(tmp_path):
    m=api();value,seal=native_positive(tmp_path,['1\tbad\tA\tC\t.\t.']*2)
    value['max_source_line_bytes']=value['source_bytes']-27
    with pytest.raises(m.Refusal):m.validate_counts(value,seal)


def test_invalid_position_raw_alt_count_cannot_exceed_source_bytes(tmp_path):
    m=api();value,seal=native_positive(tmp_path,['1\tbad\tA\tC\t.\t.'])
    value['raw_alt_tokens']=value['source_bytes']+1
    with pytest.raises(m.Refusal):m.validate_counts(value,seal)
