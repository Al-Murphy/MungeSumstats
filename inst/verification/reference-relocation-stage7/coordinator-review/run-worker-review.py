from pathlib import Path
import hashlib,json,os,subprocess,sys,time

base=Path('/private/tmp/gwas-harvest-stage2')
root=Path(__file__).resolve().parent
sys.path.insert(0,str(base/'results/verification-stage7'))
from guard_common import guarded

profile=sys.argv[1]
out=root/profile;out.mkdir(exist_ok=True);(out/'tmp').mkdir(exist_ok=True)
repo=root/'harvest'
lib=(base/'results/verification-stage7/mss-relocation-install/lib'
     if profile=='supported' else base/'Rlib/preservation-v7')
dll=lib/'MungeSumstats/libs/MungeSumstats.so'
initial=hashlib.sha256(dll.read_bytes()).hexdigest()
env=os.environ.copy()
env.update(PYTHONPATH=str(repo/'src')+os.pathsep+str(repo/'tests'),
 HARVEST_TEST_MSS_SOURCE=str(root/'mss'),HARVEST_TEST_MSS_RELOCATION=profile,
 TMPDIR=str(out/'tmp'),TMP=str(out/'tmp'),TEMP=str(out/'tmp'))
env['R_LIBS_USER']=os.pathsep.join(map(str,[lib,base/'Rlib',Path('/private/tmp/harvest-Rlib'),Path('/Users/bschilder/Library/R/arm64/4.5/library')]))
command=[str(base/'harvest/.venv/bin/python'),'-m','pytest','-q',
 str(repo/'tests/test_mss_worker.py'),str(repo/'tests/test_mss_worker_relocation.py'),
 str(root/'worker_probes.py'),f'--junitxml={out}/junit.xml']
os.chdir(repo)
result=guarded(command,out,out/'guard.json',out/'run.log',
 dict(max_tree_memory_bytes=4*1024**3,min_filesystem_free_bytes=160*1024**3,
 max_case_disk_bytes=2*1024**3,poll_seconds=.25),time.monotonic()+600,env)
assert hashlib.sha256(dll.read_bytes()).hexdigest()==initial
(out/'dll-pin.json').write_text(json.dumps({'library':str(lib),'sha256':initial,'rechecked':True},indent=2)+'\n')
print(json.dumps({'profile':profile,'status':result['status'],'dll_sha256':initial}))
