"""Probe identity aliases of pre-seal controls using owned fixtures only."""
from pathlib import Path
import json, os, sys
ROOT=Path('/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-fixes-independent-stage11')
CLONE=Path('/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-fixes-stage11')
sys.path.insert(0,str(CLONE))
import test_coordinator as t

os.umask(0o077)
root=ROOT/'hardlinked-plan-gate';root.mkdir()
m=t.api();p=t.fixture(root);expected=t.fixture_roots(root)
target=root/'GRCh37/variants.tsv';alias=root/'outside-source-plan.json';os.link(target,alias)
p['builds'][0]['source_plan']=t.spec(alias)
command=p['builds'][0]['producer_command'];command[2]=str(alias)
t.revise(p,0,'producer_guard',lambda v:v.update(command=command))
t.revise(p,1,'producer_guard',lambda v:v.update(status='running'))
opened=[];original=Path.open
def observe(path,*args,**kwargs):
    if path.samefile(target):opened.append(dict(path=str(path),mode=list(args)))
    return original(path,*args,**kwargs)
Path.open=observe
try:
    try:m.validate_seals(p,expected_publications=expected)
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    else:raise AssertionError('Malformed plan accepted')
finally:Path.open=original
value=dict(second_producer_status='running',alias_is_counted_input=os.path.samefile(alias,target),
           alias_nlink=alias.stat().st_nlink,counted_input_content_opens=opened,error=error)
(ROOT/'hardlinked-plan-gate.json').write_text(json.dumps(value,indent=2)+'\n')
print(json.dumps(value,indent=2))
