from pathlib import Path
import hashlib,json,os,subprocess,time
from guard_common import guarded
os.umask(0o077)
b=Path('/private/tmp/gwas-harvest-stage2');out=b/'results/verification-stage7/mss-relocation-install';out.mkdir();lib=out/'lib';lib.mkdir();(out/'tmp').mkdir()
repo=b/'mss-context-format3-release'
for name in ('R/preservation.R','R/projection.R','src/preservation_core.h','src/RcppExports.cpp','R/RcppExports.R'):
 actual=(repo/name).read_bytes();expected=subprocess.check_output(['git','-C',str(repo),'show','8334f1c9052ae266dc9df80a35bca07868af1423:'+name]);assert actual==expected,name
env=os.environ.copy();env['R_LIBS_USER']=':'.join(map(str,[lib,b/'Rlib',Path('/private/tmp/harvest-Rlib'),Path('/Users/bschilder/Library/R/arm64/4.5/library')]))
env.update(TMPDIR=str(out/'tmp'),TMP=str(out/'tmp'),TEMP=str(out/'tmp'))
limits={'max_tree_memory_bytes':8*1024**3,'min_filesystem_free_bytes':160*1024**3,'poll_seconds':.25,'max_case_disk_bytes':4*1024**3}
command=['/usr/bin/time','-l','R','CMD','INSTALL','--preclean','--library='+str(lib),str(repo)]
guarded(command,out,out/'install-guard.json',out/'install.log',limits,time.monotonic()+900,env)
dll=lib/'MungeSumstats/libs/MungeSumstats.so';pin={'production_revision':'8334f1c9052ae266dc9df80a35bca07868af1423','dll_sha256':hashlib.sha256(dll.read_bytes()).hexdigest(),'dll':str(dll)}
(out/'installed-pin.json').write_text(json.dumps(pin,indent=2)+'\n');print(json.dumps({'install':'passed','dll_sha256':pin['dll_sha256']}),flush=True)
subprocess.run(['python3','-B',str(b/'results/verification-stage7/run_checks-relocation.py'),'mss-relocation'],check=True)
