#!/usr/bin/env python3
"""Prepare append-only profiles; execute only with exact external root GREEN approval."""
from __future__ import annotations
import argparse,ast,csv,hashlib,io,json,os,platform,re,shutil,subprocess,sys,time
from pathlib import Path
from common import ROOT,pin,sha256,check_pin,write_once,guarded,tree_usage
from verify import Gold
STATIC=['config.json','process_memory.py','common.py','verify.py','bench.py','runner.R',
        'inspect-reference.R','selfcheck.py','README.md','approval.example.json']

def config():return json.loads((ROOT/'config.json').read_text())
def basepath(value):
 p=Path(value);return p.resolve() if p.is_absolute() else (Path(config()['base'])/p).resolve()
def source_pins(folder,revision):
 out={}
 for name in config()['production_files']:
  path=folder/name;value=pin(path)
  data=subprocess.check_output(['git','-C',str(folder),'show',revision+':'+name])
  if hashlib.sha256(data).hexdigest()!=value['sha256']:raise RuntimeError('Production checkout differs from declared revision: '+name)
  out[name]=value
 return dict(path=str(folder),production_revision=revision,observed_checkout_head=subprocess.check_output(['git','-C',str(folder),'rev-parse','HEAD'],text=True).strip(),files=out)
def library_pins(folder):
 package=folder/'MungeSumstats';files={}
 for rel in ['DESCRIPTION','NAMESPACE','R/MungeSumstats','R/MungeSumstats.rdb','R/MungeSumstats.rdx','libs/MungeSumstats.so']:
  files[rel]=pin(package/rel)
 return dict(path=str(folder),files=files,dll=files['libs/MungeSumstats.so'])
def seal():
 target=ROOT/'seal.json'
 if not target.exists():write_once(target,dict(schema=1,static_files={name:pin(ROOT/name) for name in STATIC}))
 value=json.loads(target.read_text())
 for info in value['static_files'].values():check_pin(info)
 return sha256(target)
def check_profile(profile,stage=None):
 if profile['harness_seal_sha256']!=seal():raise RuntimeError('Harness seal differs from profile')
 for group in ['mss4','mss7']:
  for info in profile['libraries'][group]['files'].values():check_pin(info)
  for info in profile['sources'][group]['files'].values():check_pin(info)
 for info in profile['reference']['files'].values():check_pin(info)
 check_pin(profile['reference']['manifest'])
 for key,s in profile['stages'].items():
  if stage is None or key==stage:check_pin(s['input'])
def preflight(args):
 cfg=config();sealsum=seal();base=Path(cfg['base']);start=time.monotonic()
 tree=tree_usage(os.getpid())
 for name in ['common.py','verify.py','bench.py','selfcheck.py','process_memory.py']:ast.parse((ROOT/name).read_text())
 for name in ['runner.R','inspect-reference.R']:
  subprocess.run(['Rscript','--vanilla','-e','invisible(parse(file=commandArgs(TRUE)[[1]]))',str(ROOT/name)],check=True,capture_output=True,text=True)
 reference=basepath(cfg['reference'])
 result=subprocess.run(['Rscript','--vanilla',str(ROOT/'inspect-reference.R'),str(reference)],check=True,capture_output=True,text=True)
 rows=list(csv.DictReader(io.StringIO(result.stdout),delimiter='\t'));meta={r['name']:r['value'] for r in rows if r['kind']=='metadata'};files={}
 for row in rows:
  if row['kind']!='file':continue
  info=pin(row['path'])
  if info['md5']!=row['md5'] or info['bytes']!=int(float(row['size'])):raise RuntimeError('Reference differs from immutable RDS receipt')
  files[row['name']]=info
 gold=Gold(files['source.tsv']['path'],files['fasta']['path'],files['fasta.fai']['path'],meta['build'])
 if len(gold.panel)!=cfg['expected_unique_panel_records']:raise RuntimeError('Panel acquisition denominator changed')
 inputs={};prior={}
 for rel in cfg['prior_receipts']:prior[rel]=pin(basepath(rel))
 for name,stage in cfg['stages'].items():
  info=pin(basepath(stage['input']))
  if info['newline_count']!=stage['rows']+1:raise RuntimeError('Pinned source physical row count changed')
  with Path(info['path']).open() as f:
   if f.readline().rstrip('\r\n')!='SNP\tCHR\tBP\tA1\tA2\tBETA\tFRQ\tSE\tP':raise RuntimeError('Input header changed')
  historical=[]
  for rel,pinned in prior.items():
   receipt=json.loads(Path(pinned['path']).read_text());metrics=receipt.get('metrics',{})
   if metrics.get('input')==info['path'] and 'input_md5' in metrics:
    if metrics['input_md5']!=info['md5']:raise RuntimeError('Input differs from previous frozen benchmark receipt')
    historical.append(rel)
  if not historical:raise RuntimeError('No prior source MD5 receipt found for '+name)
  inputs[name]=dict(stage,input=info,prior_matched_receipts=historical)
 baseline_source=basepath('mss-fasta-offsets');baseline_revision=subprocess.check_output(['git','-C',str(baseline_source),'rev-parse','HEAD'],text=True).strip()
 profile=dict(schema=1,name=args.name,created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),harness_seal_sha256=sealsum,
  preparation_only=True,execution_authorized=False,
  sources=dict(mss4=source_pins(baseline_source,baseline_revision),mss7=source_pins(basepath(args.candidate_source),args.production_revision)),
  libraries=dict(mss4=library_pins(basepath(cfg['baseline_library'])),mss7=library_pins(basepath(args.candidate_library))),
  reference=dict(directory=str(reference),build=meta['build'],version=meta['version'],format=meta['core_version'],files=files,
    source=files['source.tsv']['path'],fasta=files['fasta']['path'],fai=files['fasta.fai']['path'],manifest=pin(reference/'manifest.rds'),inventory=gold.inventory),
  stages=inputs,prior_receipts=prior,limits={k:cfg[k] for k in ['max_tree_memory_bytes','min_filesystem_free_bytes','poll_seconds']},
  preflight=dict(elapsed_seconds=time.monotonic()-start,operator_tree_measurement=tree,free_bytes=shutil.disk_usage(base).free,
    platform=platform.platform(),python=sys.version.split()[0],R=subprocess.check_output(['Rscript','--version'],text=True,stderr=subprocess.STDOUT).strip()),
  scope=cfg['workload'],cache_caveat=cfg['cache'])
 if profile['libraries']['mss4']['dll']['sha256']!=cfg['baseline_dll_sha256']:raise RuntimeError('Frozen baseline4 DLL differs')
 folder=ROOT/'profiles';folder.mkdir(exist_ok=True);path=folder/(args.name+'.json');write_once(path,profile)
 print(json.dumps(dict(profile=str(path),profile_sha256=sha256(path),preflight=profile['preflight'],reference_inventory=gold.inventory,execution_authorized=False),indent=2))

def approved(profile_path,profile,approval_path,stage):
 approval=json.loads(Path(approval_path).read_text())
 if approval.get('approved') is not True or approval.get('status')!='SCIENCE_GREEN':raise RuntimeError('Root SCIENCE_GREEN approval missing')
 if approval.get('profile_sha256')!=sha256(profile_path):raise RuntimeError('Approval does not name this immutable profile')
 if stage not in approval.get('authorized_stages',[]):raise RuntimeError('This stage has not been authorized by root')
 core=profile['sources']['mss7']['files']['src/preservation_core.h']['sha256']
 if approval.get('candidate_core_sha256')!=core:raise RuntimeError('Approval production core differs')
 review=approval['scientific_review'];review_path=Path(review['path'])
 revision=profile['sources']['mss7']['production_revision']
 expected_verdict='**Verdict: GREEN for the documented scientific contract of frozen `'+revision+'`.**'
 review_text=review_path.read_text()
 if sha256(review_path)!=review['sha256'] or not any(line.startswith(expected_verdict) for line in review_text.splitlines()) or core not in review_text:raise RuntimeError('Scientific review receipt is not an exact GREEN for the pinned source/core')
 installed=approval['installed_suite'];suite_path=Path(installed['path'])
 if sha256(suite_path)!=installed['sha256']:raise RuntimeError('Installed suite receipt changed')
 suite=json.loads(suite_path.read_text())
 if suite.get('failed_assertions')!=0 or suite.get('errors')!=0 or suite.get('passed_assertions',0)<1:raise RuntimeError('Installed suite is not passing')
 if suite.get('source_revision')!=revision:raise RuntimeError('Installed suite source revision differs from pinned production source')
 dll=suite.get('installed_dll_sha256',suite.get('dll_sha256'))
 if dll!=profile['libraries']['mss7']['dll']['sha256']:raise RuntimeError('Installed suite DLL is not the pinned candidate')
 return dict(path=str(Path(approval_path).resolve()),sha256=sha256(approval_path),review_sha256=review['sha256'],suite_sha256=installed['sha256'])

def time_metrics(log):
 text=Path(log).read_text();out={}
 if sys.platform=='darwin':
  m=re.search(r'^\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys',text,re.M)
  if m:out.update(wall_seconds=float(m[1]),user_seconds=float(m[2]),system_seconds=float(m[3]))
  for label,name in [('maximum resident set size','peak_rss_bytes'),('peak memory footprint','peak_physical_footprint_bytes')]:
   m=re.search(r'^\s*(\d+)\s+'+label,text,re.M)
   if m:out[name]=int(m[1])
 else:
  m=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)
  if m:out['peak_rss_bytes']=int(m[1])*1024
 return out

def run_stage(args):
 started=time.monotonic();profile_path=Path(args.profile).resolve();profile=json.loads(profile_path.read_text());cfg=config();stage=profile['stages'][args.stage]
 approval=approved(profile_path,profile,args.approval,args.stage)
 check_profile(profile,args.stage)
 run_parent=ROOT/'runs'/profile['name'];run_parent.mkdir(parents=True,exist_ok=True)
 if args.stage!='1M':
  previous=run_parent/'1M'/'stage-summary.json'
  if not previous.is_file() or json.loads(previous.read_text()).get('passed') is not True:raise RuntimeError('Matched1M pair must pass before10M/40M')
 out=run_parent/args.stage;out.mkdir(exist_ok=False);deadline=time.monotonic()+stage['deadline_seconds_each']*len(stage['versions']);editions={};cases={}
 try:
  for version in stage['versions']:
   case=out/version;case.mkdir();tmp=case/'tmp';tmp.mkdir();library=profile['libraries'][version];dll=library['dll'];implementation=cfg['baseline_implementation'] if version=='mss4' else cfg['candidate_implementation']
   spec=dict(version=version,library=library['path'],dll_path=dll['path'],dll_md5=dll['md5'],dll_sha256=dll['sha256'],implementation=implementation,
    production_revision=profile['sources'][version]['production_revision'],input=stage['input']['path'],reference=profile['reference']['directory'],
    output=str(case/'edition'),metrics=str(case/'metrics.json'),rows=stage['rows'],threads=cfg['threads'],batch_rows=cfg['batch_rows'],workload=cfg['workload'],cache_caveat=cfg['cache'])
   write_once(case/'spec.json',spec);limits=dict(cfg,max_case_disk_bytes=stage['max_case_disk_bytes']);env=os.environ.copy();env.update(TMPDIR=str(tmp),TMP=str(tmp),TEMP=str(tmp))
   flag='-l' if sys.platform=='darwin' else '-v';cmd=['/usr/bin/time',flag,'Rscript','--vanilla',str(ROOT/'runner.R'),str(case/'spec.json')]
   receipt=guarded(cmd,case,case/'process-guard.json',case/'process.log',limits,min(deadline,time.monotonic()+stage['deadline_seconds_each']),env)
   metrics=json.loads((case/'metrics.json').read_text());external=time_metrics(case/'process.log')
   if external.get('peak_rss_bytes',0)>cfg['max_tree_memory_bytes'] or external.get('peak_physical_footprint_bytes',0)>cfg['max_tree_memory_bytes']:raise RuntimeError('External time peak exceeds8GiB target')
   cases[version]=dict(metrics=metrics,external_time=external,guard=receipt);editions[version]=str(case/'edition')
  command=[sys.executable,'-B',str(ROOT/'verify.py'),'--profile',str(profile_path),'--stage',args.stage,'--receipt',str(out/'correctness.json')]
  for version,folder in editions.items():command.extend(['--edition',version+'='+folder])
  limits=dict(cfg,max_case_disk_bytes=stage['max_case_disk_bytes']*len(stage['versions']))
  verification=guarded(command,out,out/'verification-guard.json',out/'verification.log',limits,deadline)
  check_profile(profile,args.stage)
  correctness=json.loads((out/'correctness.json').read_text())
  result=dict(schema=1,stage=args.stage,passed=correctness['passed'],profile_sha256=sha256(profile_path),approval=approval,cases=cases,
   correctness=correctness,verification_guard=verification,total_operator_wall_seconds=time.monotonic()-started,
   processing_reference_verify_package_setup_external_wall_are_separate=True,workload=cfg['workload'],cache_caveat=cfg['cache'],not_unique_variant_count=True)
  write_once(out/'stage-summary.json',result);print(json.dumps(result,indent=2))
 except BaseException as exc:
  write_once(out/'failed-stage.json',dict(stage=args.stage,passed=False,error=type(exc).__name__+': '+str(exc),profile_sha256=sha256(profile_path),elapsed_seconds=time.monotonic()-started,natural_oom=False))
  raise

def main():
 ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True);cfg=config()
 p=sub.add_parser('preflight');p.add_argument('--name',required=True);p.add_argument('--candidate-library',default=cfg['candidate_library_default']);p.add_argument('--candidate-source',default=cfg['candidate_source_default']);p.add_argument('--production-revision',default=cfg['candidate_production_revision_default'])
 p=sub.add_parser('run-stage');p.add_argument('--profile',required=True);p.add_argument('--approval',required=True);p.add_argument('--stage',choices=list(cfg['stages']),required=True)
 a=ap.parse_args()
 if a.action=='preflight':
  if not re.fullmatch('[a-zA-Z0-9_-]+',a.name):raise RuntimeError('Unsafe profile name')
  preflight(a)
 else:run_stage(a)
if __name__=='__main__':main()
