"""Observe guard ordering and preserve legitimate hardlinked-data compatibility."""
from pathlib import Path
import json, os, sys
ROOT=Path('/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-fixes-independent-stage11')
CLONE=Path('/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-fixes-stage11')
sys.path.insert(0,str(CLONE))
import test_coordinator as t

os.umask(0o077);m=t.api();results={}
root=ROOT/'guard-first-order';root.mkdir();p=t.fixture(root)
t.revise(p,1,'producer_guard',lambda v:v.update(status='running'))
roles={b[role]['path']:b['build']+' '+role for b in p['builds'] for role in ['producer_guard','source_plan','receipt']}
order=[];original=Path.open
def observe(path,*args,**kwargs):
    if str(path) in roles:order.append(roles[str(path)])
    return original(path,*args,**kwargs)
Path.open=observe
try:
    try:m.validate_seals(p,expected_publications=t.fixture_roots(root))
    except m.Refusal as exc:error=str(exc)
    else:raise AssertionError('Running second guard accepted')
finally:Path.open=original
results['running_second_guard_control_open_order']={'order':order,'refusal':error}

root=ROOT/'hardlinked-data-positive';root.mkdir();p=t.fixture(root)
target=root/'GRCh37/variants.tsv';alias=root/'legitimate-data-copy-link.tsv';os.link(target,alias)
before={str(x):t.spec(x) for b in p['builds'] for x in Path(b['final_dir']).iterdir()}
result=m.count_worker(p,root/'case',expected_publications=t.fixture_roots(root))
after={str(x):t.spec(x) for b in p['builds'] for x in Path(b['final_dir']).iterdir()}
assert before==after
results['hardlinked_counted_data_positive']={'status':result['status'],'data_nlink':target.stat().st_nlink,
                  'publication_files_unchanged':before==after,'full_reference_admitted':result['full_reference_admitted'],
                  'index_builds_executed':result['index_builds_executed']}
(ROOT/'guard-order-and-data-link.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
