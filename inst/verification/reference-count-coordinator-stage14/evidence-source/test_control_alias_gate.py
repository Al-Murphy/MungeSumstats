"""Owned tiny fixtures for the remaining small-control identity/order gate."""
import os
from pathlib import Path

import pytest
from test_coordinator import api, fixture, fixture_roots, revise, spec


@pytest.mark.parametrize('role',['source_plan','producer_guard'])
@pytest.mark.parametrize('build_index',[0,1])
def test_small_control_alias_never_opens_counted_identity(tmp_path,monkeypatch,role,build_index):
    m=api();p=fixture(tmp_path)
    revise(p,1,'producer_guard',lambda v:v.update(status='running'))
    target=Path(p['builds'][build_index]['final_dir'])/'variants.tsv'
    alias=tmp_path/('outside-'+role+'.json');os.link(target,alias)
    p['builds'][build_index][role]=spec(alias)
    if role=='source_plan':
        command=p['builds'][build_index]['producer_command'];command[2]=str(alias)
        revise(p,build_index,'producer_guard',lambda v:v.update(command=command))
    opened=[];original=Path.open
    def observe(path,*args,**kwargs):
        if path.samefile(target):opened.append(str(path))
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
    assert opened==[]


@pytest.mark.parametrize('failure',['running','wrong_command'])
def test_both_guards_authenticate_before_plan_receipt_or_artifact_content(tmp_path,monkeypatch,failure):
    m=api();p=fixture(tmp_path)
    revise(p,1,'producer_guard',lambda v:v.update(status='running') if failure=='running' else v.update(command=['wrong']))
    forbidden={b[k]['path'] for b in p['builds'] for k in ('source_plan','receipt')}
    forbidden.update(s['path'] for s in p['artifacts'].values())
    opened=[];original=Path.open
    def observe(path,*args,**kwargs):
        opened.append(str(path));assert str(path) not in forbidden,'content before both guard seals: '+str(path)
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
    assert set(opened)=={b['producer_guard']['path'] for b in p['builds']}


def test_both_successful_guards_precede_all_other_control_content(tmp_path,monkeypatch):
    m=api();p=fixture(tmp_path)
    guards={b['producer_guard']['path'] for b in p['builds']};opened=[];original=Path.open
    def observe(path,*args,**kwargs):opened.append(str(path));return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    seals=m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
    assert [s['build'] for s in seals]==['GRCh37','GRCh38']
    assert max(i for i,path in enumerate(opened) if path in guards)<min(i for i,path in enumerate(opened) if path not in guards)
    assert guards.issubset(opened)


@pytest.mark.parametrize('name',['variants.tsv','reference.fa','contig_aliases.tsv'])
def test_hardlinked_counted_data_remains_valid_without_copy(tmp_path,name):
    m=api();p=fixture(tmp_path);target=tmp_path/'GRCh37'/name
    alias=tmp_path/('data-link-'+name);os.link(target,alias)
    before={str(f):spec(f) for b in p['builds'] for f in Path(b['final_dir']).iterdir()}
    result=m.count_worker(p,tmp_path/'case',expected_publications=fixture_roots(tmp_path))
    assert before=={str(f):spec(f) for b in p['builds'] for f in Path(b['final_dir']).iterdir()}
    assert target.samefile(alias) and target.stat().st_nlink==2
    assert [b['counts']['source_rows'] for b in result['builds']]==[1,1]
    assert result['full_reference_admitted'] is False and result['index_builds_executed']==0


@pytest.mark.parametrize('role',['source_plan','producer_guard'])
def test_single_link_rule_is_limited_to_small_external_controls(tmp_path,monkeypatch,role):
    m=api();p=fixture(tmp_path);control=Path(p['builds'][0][role]['path'])
    os.link(control,tmp_path/('extra-'+role+'.json'))
    original=Path.open;opened=[]
    def observe(path,*args,**kwargs):opened.append(str(path));return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',observe)
    with pytest.raises(m.Refusal):m.validate_seals(p,expected_publications=fixture_roots(tmp_path))
    assert opened==[]
