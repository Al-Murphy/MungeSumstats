"""Private bounded final integrations; no provider deployment or billing action."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,urllib.parse
from guard_common import guarded

os.umask(0o077)
base=Path('/private/tmp/gwas-harvest-stage2');phase=sys.argv[1]
if phase not in ('mss-relocation',):raise ValueError('Explicit verification phase required')
root=Path(__file__).resolve().parent;out=root/phase;out.mkdir();(out/'tmp').mkdir()
env=os.environ.copy();env.update(TMPDIR=str(out/'tmp'),TMP=str(out/'tmp'),TEMP=str(out/'tmp'))
private={}
if phase=='harvest-mss7':
 repo=base/'harvest';private=json.loads((base/'results/metadata-pg/test.env.json').read_text());env.update(private)
 env.update(PYTHONPATH=str(repo/'src'),HARVEST_TEST_MSS_SOURCE=str(base/'mss-orientation-v7'),MSS_TEST_RLIB=str(base/'Rlib'))
 libraries=[base/'Rlib/preservation-v7',base/'Rlib/reference-contigs-v1',base/'Rlib',Path('/private/tmp/harvest-Rlib'),Path('/Users/bschilder/Library/R/arm64/4.5/library')]
 command=['/usr/bin/time','-l',str(repo/'.venv/bin/python'),'-m','pytest','-q']
else:
 repo=base/'mss-context-format3-release';libraries=[base/'results/verification-stage7/mss-relocation-install/lib',base/'Rlib',Path('/private/tmp/harvest-Rlib'),Path('/Users/bschilder/Library/R/arm64/4.5/library')]
 command=['/usr/bin/time','-l','Rscript','--vanilla',str(root/'full_mss.R'),str(repo),str(out),str(libraries[0]/'MungeSumstats')]
env['R_LIBS_USER']=':'.join(map(str,libraries))
files=subprocess.check_output(['git','-C',str(repo),'ls-files','src','R','tests','migrations'],text=True).splitlines()
files=[name for name in files if (repo/name).is_file()]
pins={name:hashlib.sha256((repo/name).read_bytes()).hexdigest() for name in files}
state={'phase':phase,'checkout_head':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),'source_pins':pins,'libraries':list(map(str,libraries))}
if phase=='mss-relocation':
 state['production_revision']='8334f1c9052ae266dc9df80a35bca07868af1423'
 dll=libraries[0]/'MungeSumstats/libs/MungeSumstats.so';state['dll_sha256']=hashlib.sha256(dll.read_bytes()).hexdigest()
 if state['dll_sha256']!=json.loads((base/'results/verification-stage7/mss-relocation-install/installed-pin.json').read_text())['dll_sha256']:raise ValueError('Relocation installed DLL differs')
else:
 dll=libraries[0]/'MungeSumstats/libs/MungeSumstats.so';state['dll_sha256']=hashlib.sha256(dll.read_bytes()).hexdigest()
 if state['dll_sha256']!='c1c0a8b29479c33d8b4cfa4646e46e39c4da436686f4ca2da12a402251003d50':raise ValueError('Scientific7 DLL differs')
(out/'source-pins.json').write_text(json.dumps(state,indent=2)+'\n')
limits={'max_tree_memory_bytes':8*1024**3,'min_filesystem_free_bytes':160*1024**3,'poll_seconds':.25,'max_case_disk_bytes':4*1024**3}
old=Path.cwd();os.chdir(repo)
try:
 receipt=guarded(command,out,out/'guard.json',out/'run.private.log',limits,time.monotonic()+900,env)
finally:
 os.chdir(old)
 text=(out/'run.private.log').read_text()
 secrets=[]
 for value in private.values():
  if isinstance(value,str) and value:
   secrets.append(value)
   try:
    parsed=urllib.parse.urlsplit(value)
    if parsed.password:secrets.append(urllib.parse.unquote(parsed.password))
   except ValueError:pass
 for value in sorted(secrets,key=len,reverse=True):text=text.replace(value,'<private fixture configuration>')
 (out/'run.log').write_text(text)
for name,digest in pins.items():
 if hashlib.sha256((repo/name).read_bytes()).hexdigest()!=digest:raise ValueError('Source changed during verification: '+name)
if hashlib.sha256(dll.read_bytes()).hexdigest()!=state['dll_sha256']:raise ValueError('Installed DLL changed during verification')
print(json.dumps({'phase':phase,'status':receipt['status'],'source_files_rechecked':len(pins),'public_log':str(out/'run.log')}),flush=True)
