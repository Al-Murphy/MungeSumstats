"""Actual frozen v4/v5 tiny publications and existing count-only binary; no writer."""
import json
from pathlib import Path

import pytest
from test_coordinator import spec, fixture_roots
from test_coordinator_integration import produce_both, run_profile


@pytest.mark.parametrize('versions,schema,empty', [
    (('4','4'), 1, False), (('5','5'), 2, True), (('4','5'), 2, True),
])
def test_actual_frozen_preparers_and_guarded_count(tmp_path, versions, schema, empty):
    p=produce_both(tmp_path,versions,schema_version=schema,empty_history=empty)
    before={str(x):spec(x) for b in p['builds'] for x in Path(b['final_dir']).iterdir()}
    result=run_profile(p,tmp_path,'guarded-count',expected_publications=fixture_roots(tmp_path,published=True))
    assert result['guard']['status']=='success'
    assert result['decision']=='no_admission' and result['full_reference_admitted'] is False
    assert result['index_builds_executed']==0
    for b,r,version in zip(p['builds'],result['builds'],versions,strict=True):
        receipt=json.loads(Path(b['receipt']['path']).read_text())
        assert receipt['preparation_version']=='biodb-reference-stream-'+version
        assert r['counts']['source_rows']==b['expected_emitted_rows']==(4 if empty and version=='5' else 3)
        text=(Path(b['final_dir'])/'variants.tsv').read_text()
        assert 'C,G,AA,<DEL>,*' in text and 'AA,AAA' in text
        if schema==2:
            assert r['preparation_version']==b['preparation_version']
            assert r['producer_source']['sha256']==b['producer_source']['sha256']
        if empty and version=='5':
            assert receipt['alias_history_status']=='complete_supplied_history'
            assert receipt['counts']['merged_empty_records']==1
    assert before=={str(x):spec(x) for b in p['builds'] for x in Path(b['final_dir']).iterdir()}
    assert not any((Path(b['final_dir'])/'fasta.fai').exists() for b in p['builds'])
