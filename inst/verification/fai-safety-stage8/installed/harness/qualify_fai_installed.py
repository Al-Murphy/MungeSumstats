from pathlib import Path
import hashlib,io,json,os,subprocess,tarfile,time
from guard_common import guarded

os.umask(0o077)
base=Path('/private/tmp/gwas-harvest-stage2')
revision='03d11e92b417add0183dea4fb77716f5502d5735'
out=base/'results/verification-stage8/mss-fai-installed';out.mkdir(parents=True)
source=out/'source';source.mkdir();lib=out/'lib';lib.mkdir();tmp=out/'tmp';tmp.mkdir()
payload=subprocess.check_output(['git','-C',str(base/'mss-fai-safety-stage8'),'archive',revision])
with tarfile.open(fileobj=io.BytesIO(payload)) as archive:archive.extractall(source,filter='data')
names=subprocess.check_output(['git','-C',str(base/'mss-fai-safety-stage8'),'ls-tree','-r','--name-only',revision],text=True).splitlines()
pins={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in names if (source/name).is_file()}
assert pins['src/preservation_core.h']=='4affcc7bdb761b1c86d22e2f72eb079434fbb2404d7f5c4d257fda2de52c6ddc'
env=os.environ.copy();env.update(TMPDIR=str(tmp),TMP=str(tmp),TEMP=str(tmp))
env['R_LIBS_USER']=':'.join(map(str,[lib,base/'Rlib',Path('/private/tmp/harvest-Rlib'),Path('/Users/bschilder/Library/R/arm64/4.5/library')]))
limits={'max_tree_memory_bytes':8*1024**3,'min_filesystem_free_bytes':160*1024**3,
        'max_case_disk_bytes':4*1024**3,'poll_seconds':.25}
guarded(['/usr/bin/time','-l','R','CMD','INSTALL','--preclean','--library='+str(lib),str(source)],
        out,out/'install-guard.json',out/'install.log',limits,time.monotonic()+900,env)
dll=lib/'MungeSumstats/libs/MungeSumstats.so'
pin={'revision':revision,'source_pins':pins,'dll_sha256':hashlib.sha256(dll.read_bytes()).hexdigest(),
     'source_archive_sha256':hashlib.sha256(payload).hexdigest(),'library':str(lib)}
(out/'installed-pin.json').write_text(json.dumps(pin,indent=2)+'\n')
suite=out/'suite';suite.mkdir();(suite/'tmp').mkdir()
env.update(TMPDIR=str(suite/'tmp'),TMP=str(suite/'tmp'),TEMP=str(suite/'tmp'))
command=['/usr/bin/time','-l','Rscript','--vanilla',str(base/'results/verification-stage7/full_mss.R'),str(source),str(suite),str(lib/'MungeSumstats')]
guarded(command,suite,suite/'guard.json',suite/'run.log',limits,time.monotonic()+900,env)
assert all(hashlib.sha256((source/name).read_bytes()).hexdigest()==sha for name,sha in pins.items())
assert hashlib.sha256(dll.read_bytes()).hexdigest()==pin['dll_sha256']
print(json.dumps({'status':'passed','revision':revision,'dll_sha256':pin['dll_sha256'],'source_files_rechecked':len(pins)}))
