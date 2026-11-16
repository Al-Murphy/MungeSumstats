import json,os,subprocess
from pathlib import Path
b=Path('/private/tmp/gwas-harvest-stage2');source=b/'mss-orientation-v7';out=b/'results/mss-preservation-v7/full-suite'
env=os.environ.copy();env['R_LIBS_USER']=':'.join(map(str,[b/'Rlib/preservation-v7',b/'Rlib/fasta-offsets-v1',b/'Rlib/reference-contigs-v1',b/'Rlib',Path('/private/tmp/harvest-Rlib'),Path('/Users/bschilder/Library/R/arm64/4.5/library')]))
commands=[(['R','CMD','INSTALL','--preclean','--library='+str(b/'Rlib/preservation-v7'),str(source)],'install.log'),(['python3',str(b/'results/full-reference/guard.py'),'--receipt',str(out/'guard.json'),'--filesystem',str(b),'--','/usr/bin/time','-l','Rscript','--vanilla',str(out/'run-full-suite.R')],'full-suite.log')]
for command,name in commands:
 with (out/name).open('x') as log:
  result=subprocess.run(command,cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT)
 print(json.dumps({'phase':name,'exit_code':result.returncode}),flush=True)
 if result.returncode:raise SystemExit(result.returncode)
