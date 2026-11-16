"""Independent review probes, exclusively new owned tiny fixtures."""
from pathlib import Path
import copy, hashlib, importlib.util, json, os, subprocess, sys

ROOT=Path('/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-independent-stage11')
CLONE=Path('/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-stage11')
sys.path.insert(0,str(CLONE))
import test_coordinator as t
import test_coordinator_integration as it

def own(name):
    p=ROOT/name;p.mkdir();return p

def main():
    os.umask(0o077)
    findings={}
    m=t.api()
    # Ordinary fixture positive control, exact same trusted source/count binary.
    root=own('positive');p=t.fixture(root)
    r=m.count_worker(p,root/'case')
    findings['ordinary_positive']={'status':r['status'],'counts':[x['counts'] for x in r['builds']],
                                  'full_reference_admitted':r['full_reference_admitted']}
    # The malformed response is a native-output fault injection, following the
    # author's own test strategy, not a claim to sandbox a malicious operator.
    root=own('impossible-counts');p=t.fixture(root)
    fake=root/'corrupt-native'
    fake.write_text('#!'+sys.executable+'\nimport json,subprocess,sys\n'+
                    'v=json.loads(subprocess.check_output(['+repr(str(CLONE/'count_reference'))+',*sys.argv[1:]]))\n'+
                    'for k in ["records_bytes","keys_bytes","raw_sort_bytes","raw_key_rows","final_key_rows","offsets_bytes","rejected_bytes","fai_upper_bytes","raw_alt_tokens"]:v[k]=0\n'+
                    'print(json.dumps(v))\n')
    fake.chmod(0o700);p['artifacts']['native']=t.spec(fake)
    r=m.count_worker(p,root/'case')
    findings['impossible_native_counts_accepted']={'status':r['status'],
         'builds':r['builds'],'worker_receipt_written':(root/'case'/'worker.json').exists()}
    # Count artifacts must stay outside immutable publications.
    root=own('case-in-publication');p=t.fixture(root)
    final=Path(p['builds'][0]['final_dir']);case=final/'new-count-evidence'
    before=sorted(x.name for x in final.iterdir())
    r=m.count_worker(p,case)
    findings['case_inside_publication_accepted']={'status':r['status'],
          'publication_before':before,'publication_after':sorted(x.name for x in final.iterdir()),
          'case_files':sorted(x.name for x in case.iterdir())}
    # Correctly formed trusted profile stored inside publication is accepted.
    root=own('profile-in-publication');p=t.fixture(root)
    profile=Path(p['builds'][0]['final_dir'])/'operator-profile.json'
    pp=t.write(profile,p)
    loaded,_,_=m.load_profile(profile,pp['sha256'])
    findings['profile_inside_publication_accepted']={'loaded_context':loaded['context'],
          'profile':str(profile)}
    # Reject artifact control in a publication before opening any input.
    root=own('gate-mutation');p=t.fixture(root)
    target=Path(p['builds'][0]['final_dir'])/'variants.tsv'
    p['artifacts']['native']=t.spec(target)
    t.revise(p,1,'producer_guard',lambda v:v.update(status='running'))
    opened=[];original=Path.open
    def observe(path,*args,**kwargs):
        if path==target:opened.append(str(path))
        return original(path,*args,**kwargs)
    Path.open=observe
    try:
        try:m.validate_seals(p)
        except m.Refusal as exc:reason=str(exc)
        else:raise AssertionError('Invalid artifact gate accepted')
    finally:Path.open=original
    findings['gate_mutation_refused_without_input_open']={'reason':reason,'input_opens':opened}
    # Independently generated randomized exact arithmetic and one-byte boundaries.
    import random
    rng=random.Random(20261116)
    for i in range(300):
        def one():
            ts,rs,us,ks,ns,js,fs,als,man=[rng.randrange(10**6) for _ in range(9)]
            c=fs+als+man+19
            expected={'small_files_upper_bytes':c,'retained_bytes':ts+rs+ks+8*ns+js+c,
                      'phase_bytes':max(ts+rs+ks+8*ns+js+c,rs+2*us+js+c,rs+us+ks+8*ns+js+c)}
            got=m.logical_size(t.counts(source_bytes=ts,records_bytes=rs,raw_sort_bytes=us,
                  keys_bytes=ks,final_key_rows=ns,rejected_bytes=js,fai_upper_bytes=fs),als,man)
            assert got==expected
            return got
        a,b=one(),one();allocation=rng.randrange(1000);concurrent=rng.randrange(1000)
        peaks=[max(a['phase_bytes'],a['retained_bytes']+b['phase_bytes'],a['retained_bytes']+b['retained_bytes']),
               max(b['phase_bytes'],b['retained_bytes']+a['phase_bytes'],a['retained_bytes']+b['retained_bytes'])]
        for index,order in enumerate(['GRCh37_then_GRCh38','GRCh38_then_GRCh37']):
            boundary=peaks[index]+170*1024**3+allocation+concurrent
            assert m.capacity([a,b],boundary,allocation,concurrent)['orders_fit_with_reserves'][order]
            assert not m.capacity([a,b],boundary-1,allocation,concurrent)['orders_fit_with_reserves'][order]
    findings['arithmetic']={'random_pairs':300,'order_boundary_checks':1200,'status':'passed'}
    (ROOT/'independent-probes.json').write_text(json.dumps(findings,indent=2)+'\n')
    print(json.dumps({k:(v.get('status','observed')) for k,v in findings.items()},indent=2))

if __name__=='__main__':main()
