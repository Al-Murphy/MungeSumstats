"""Bounded FAI scientific/portable qualification; never a throughput run."""
from pathlib import Path
import hashlib,json,re,shutil,subprocess

BASE=Path('/private/tmp/gwas-harvest-stage2')
ROOT=BASE/'results/mss-fai-safety-stage8'
SOURCE=BASE/'mss-fai-safety-stage8'
OLD=BASE/'review-mss-context-format3-stage7'
WORK=ROOT/'qualification';WORK.mkdir()
def run(command,log):
 with log.open('wb') as stream:
  p=subprocess.run(list(map(str,command)),stdout=stream,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError(str(command)+'; exit '+str(p.returncode)+'; '+str(log))
def sha(path):
 with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

driver=WORK/'driver.cpp';driver.write_bytes((OLD/'driver.cpp').read_bytes())
run(['clang++','-std=c++17','-O2','-I'+str(SOURCE/'src'),driver,'-o',WORK/'core'],WORK/'compile.log')
summary={}
for fmt in [2,3]:
 d=WORK/f'format{fmt}';d.mkdir();(d/'core').symlink_to(WORK/'core')
 for name in ['oracle.py','haplotype-contract-oracle.py','contract-gold.py']:
  text=(OLD/f'format{fmt}'/name).read_text().replace('ece4fcdfad4aa9c2ce85434f67694190872442dd','6d6ec4187c69e113eda28fd2d68d986254a330c4')
  if name=='contract-gold.py':
   text=text.replace(str(OLD/f'format{fmt}/old-fixtures'),str(d/'old-fixtures'))
  (d/name).write_text(text)
 fixtures=d/'old-fixtures';fixtures.mkdir()
 for case in ['unknown-anchor','boundary','alternate-background-repeat','singlebase-alternate-indel-background']:
  out=fixtures/case;out.mkdir();old=OLD/f'format{fmt}/old-fixtures'/case
  for name in ['tiny.fa','sites.tsv']:(out/name).write_bytes((old/name).read_bytes())
  run([WORK/'core',f'prepare{fmt}',out/'tiny.fa',out/'sites.tsv',out/'index'],out/'prepare.log')
 for name in ['oracle.py','haplotype-contract-oracle.py','contract-gold.py']:
  run(['python3','-B',d/name],d/(name+'.log'))
 records={}
 for name in ['plain-oracle-output.tsv','left-oracle-output.tsv','right-oracle-output.tsv','both-oracle-output.tsv','haplotype-output.tsv']:
  assert (d/name).read_bytes()==(OLD/f'format{fmt}'/name).read_bytes(),(fmt,name)
  records[name]=dict(rows=(d/name).read_bytes().count(b'\n'),sha256=sha(d/name))
 summary[f'format{fmt}']=dict(oracle=json.loads((d/'oracle-summary.json').read_text()),haplotype=json.loads((d/'haplotype-summary.json').read_text()),explicit=json.loads((d/'contract-gold-summary.json').read_text()),complete_native_output_exact_to_science7=records)
for runner in ['preservation.R','registered.R']:
 env={'MSS_TEST_RLIB':'/private/tmp/harvest-Rlib'}
 import os
 fullenv=os.environ.copy();fullenv.update(env)
 with (WORK/(runner+'.log')).open('wb') as log:
  p=subprocess.run(['Rscript',str(SOURCE/'tests/standalone'/runner),str(SOURCE)],env=fullenv,stdout=log,stderr=subprocess.STDOUT)
 assert p.returncode==0,runner
 text=(WORK/(runner+'.log')).read_text();counts=[int(x) for x in re.findall(r'Test passed with (\d+) successes',text)];singles=len(re.findall(r'Test passed with 1 success',text))
 assert not re.search(r'Test failed|Failure|Error|Warning',text),runner
 summary[runner]=dict(blocks=len(counts)+singles,assertions=sum(counts)+singles,failures=0,errors=0,warnings=0)
(ROOT/'qualification-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v if k.endswith('.R') else {'SNP_rows':v['oracle'].get('queries',v['oracle'].get('counts')), 'explicit':v['explicit']['counts']} for k,v in summary.items()}))
