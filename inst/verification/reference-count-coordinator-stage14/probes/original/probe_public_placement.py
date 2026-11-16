"""Public launch path placement probes on new owned fixtures only."""
from pathlib import Path
import json, os, sys
ROOT=Path('/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-independent-stage11')
CLONE=Path('/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11')
sys.path.insert(0,str(CLONE))
import test_coordinator as t
import test_coordinator_integration as it

os.umask(0o077)
suffix=sys.argv[1] if len(sys.argv)>1 else ''
m=t.api();observations={}
for name,inside_profile,inside_case in [('public-positive',False,False),
                                       ('public-profile-in-publication',True,False),
                                       ('public-case-in-publication',False,True)]:
    root=ROOT/(name+suffix);root.mkdir()
    p=it.actual_profile(root)
    final=Path(p['builds'][0]['final_dir'])
    profile=(final if inside_profile else root)/'operator-profile.json'
    pp=t.write(profile,p)
    case=(final if inside_case else root)/'count-case'
    before=sorted(x.name for x in final.iterdir())
    try:
        result=m.launch(profile,pp['sha256'],case)
        status=result['status'];error=None
    except Exception as exc:
        status='exception';error=type(exc).__name__+': '+str(exc)
    observations[name]={'status':status,'error':error,
       'publication_before':before,'publication_after':sorted(x.name for x in final.iterdir()),
       'case_exists':case.exists(),'summary_exists':(case/'summary.json').exists()}

root=ROOT/('oversize-top-level-control'+suffix);root.mkdir();target=root/'oversize.json'
with target.open('xb') as handle:handle.truncate(16*1024**2+1)
pp=t.spec(target);original=Path.open;opens=[]
def observe(path,*args,**kwargs):
    if path==target:opens.append({'args':list(args)})
    return original(path,*args,**kwargs)
Path.open=observe
try:
    try:m.load_profile(target,pp['sha256'])
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    else:raise AssertionError('Oversize control accepted')
finally:Path.open=original
observations['oversize_control_hashed_before_size_refusal']={'bytes':target.stat().st_size,
                      'opens_after_starting_load_profile':opens,'error':error}
(ROOT/('public-placement-probes'+suffix+'.json')).write_text(json.dumps(observations,indent=2)+'\n')
print(json.dumps(observations,indent=2))
