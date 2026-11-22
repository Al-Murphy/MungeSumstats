"""Independent local reconciliation of the completed comparison; no workload rerun."""
import csv,hashlib,json,math,statistics
from decimal import Decimal
from pathlib import Path
HERE=Path('/private/tmp/mss-combined-results-20261122');OP=Path('/private/tmp/mss-combined-operator-20261122');RET=OP/'returned'
read=lambda p:json.loads(p.read_bytes())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def tsv(p):
 with p.open(newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))
def close(a,b):
 if a is None or b is None:assert a is b
 else:assert math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-8),(a,b)
copy=read(OP/'controller-copy-proof.json');assert len(copy['controller_pins'])==88 and not copy['absent']
for n,v in copy['controller_pins'].items():assert (RET/n).stat().st_size==v['bytes'] and sha(RET/n)==v['sha256']
for root,index,count in [(OP,OP/'SOURCE-PINS.json',24),(Path('/private/tmp/mss-combined-installed-source-20261122'),OP/'PAYLOAD-PINS.json',17)]:
 pins=read(index);assert len(pins)==count
 for n,v in pins.items():assert (root/n).stat().st_size==v['bytes'] and sha(root/n)==v['sha256']
a=read(HERE/'analysis.json');r=read(RET/'output/comparison-receipt.json');goldpins=read(Path('/private/tmp/mss-combined-installed-source-20261122/retained-pins.json'))['golds']
cases=tsv(RET/'output/functional/receipt/tests.tsv');q=read(RET/'output/functional/verified.json')
assert len(cases)==q['tests']==269 and sum(int(c['passed']) for c in cases)==q['passed']==3077
assert all(c['skipped']==c['error']=='FALSE' and int(c['failed'])==0 for c in cases)
assert sum(int(c['warning']) for c in cases)==2 and tsv(RET/'output/functional/receipt/conditions.tsv')==q['conditions']
assert [(c['file'],c['test'],int(c['passed'])) for c in cases]==[(c['file'],c['test'],c['passed']) for c in q['cases']]
assert q==r['qualified_candidate']==read(Path('/private/tmp/mss-file-batch-lifetime-operator-20261122/returned/output/functional/verified.json'))
scalar={'operation_seconds':'processing.elapsed','operation_plus_closing_gc_seconds':'processing_plus_closing_gc.elapsed','closing_gc_seconds':'closing_gc.elapsed','operation_gc_elapsed_seconds':'gc_operation_delta.3','operation_plus_closing_gc_collection_elapsed_seconds':'gc_operation_plus_closing_delta.3','gc_through_pre_receipt_elapsed_seconds':'gc_through_receipt_delta.3'}
gc_cpu={'operation_gc_cpu_seconds':'gc_operation_delta','operation_plus_closing_gc_cpu_seconds':'gc_operation_plus_closing_delta','gc_through_pre_receipt_cpu_seconds':'gc_through_receipt_delta'}
def extract(key,m,p):
 if key in scalar:return Decimal(m[scalar[key]])
 if key=='worker_wall_seconds':return Decimal(str(p['wall_seconds']))
 if key=='worker_cpu_seconds':return Decimal(str(p['user_seconds']))+Decimal(str(p['system_seconds']))
 if key=='peak_rss_bytes':return Decimal(p['peak_rss_bytes'])
 if key in gc_cpu:return sum(Decimal(m[gc_cpu[key]+'.'+str(i)]) for i in [1,2,4,5])
 phase,field=key.rsplit('.',1);i=['user_cpu','system_cpu','elapsed','child_user_cpu','child_system_cpu'].index(field)+1
 return Decimal(m[phase+'.'+str(i)])
review={};comparisons=0
for label,rows in [('1M',1000000),('10M',10000000)]:
 records={};rosters=[]
 for rec in r['results'][label]['measurements']:
  folder=RET/'output/runs'/label/rec['run'];m={x['key']:x['value'] for x in tsv(folder/'metrics.tsv')};p=read(folder/'process.json')
  assert m==rec['metrics'] and p==rec['process'] and p['exit_code']==0 and m['threads']=='1'
  assert float(m['source_rows'])==float(m['normalized_rows'])==rows
  assert rec['role']==m['role']==rec['run'][-1];records[rec['run']]=(m,p)
  rosters.append([{k:v for k,v in x.items() if k!='path'} if x['package']=='MungeSumstats' else x for x in tsv(folder/'namespaces.tsv')])
 assert len(records)==6 and all(x==rosters[0] for x in rosters)
 prior=Path('/private/tmp/mss-file-batch-lifetime-operator-20261122/returned/output/runs/pair-01-B/namespaces.tsv') if label=='1M' else Path('/private/tmp/mss-file-batch-lifetime-10M-operator-20261122/returned/output/runs/pair-01-B/namespaces.tsv')
 prior_roster=[{k:v for k,v in x.items() if k!='path'} if x['package']=='MungeSumstats' else x for x in tsv(prior)]
 assert prior_roster==rosters[0]
 summary=a['results'][label];assert len(summary['medians'])==47 and len(summary['pairs'])==3
 all_values={key:[] for key in summary['medians']}
 for pair in [1,2,3]:
  A=records[f'pair-{pair:02d}-A'];B=records[f'pair-{pair:02d}-B'];expected=summary['pairs'][pair-1]
  assert expected['pair']==pair and expected['order']==['AB','BA','AB'][pair-1]
  for key in summary['medians']:
   av=extract(key,*A);bv=extract(key,*B);dv=bv-av;pc=Decimal(100)*dv/av if av else None
   observed=expected['values'][key]
   for k,v in [('baseline',av),('candidate',bv),('difference',dv),('percent',pc)]:close(v,observed[k])
   all_values[key].append((av,bv,dv,pc));comparisons+=1
 for key,vals in all_values.items():
  for target,idx in [('baseline',0),('candidate',1),('paired_difference',2),('paired_percent',3)]:
   selected=[v[idx] for v in vals if v[idx] is not None];close(statistics.median(selected) if selected else None,summary['medians'][key][target])
 digests=read(RET/f'output/verification/{label}/artifact-digests.json');assert len(digests)==7 and str(goldpins[label]['run']) in digests
 for name,gold in goldpins[label]['edition_files'].items():
  if name!='manifest.rds':assert all(d[name]==gold for d in digests.values())
 assert (RET/f'output/verification/{label}/stable-manifests.txt').read_bytes()==b'pass\n'
 oracle=read(RET/f'output/verification/{label}/whole-SNP-oracle.json')
 assert oracle['passed'] is True and oracle['all_rows_examined'] is True and oracle['unaccounted_deltas']==0
 assert oracle['rows_expected']==oracle['versions']['mss7']['counts']['all_19_fields_correct_rows']==rows
 review[label]={k:float(summary['medians'][k]['paired_percent']) for k in ['operation_plus_closing_gc_seconds','worker_cpu_seconds','worker_wall_seconds','peak_rss_bytes']}
result={'status':'independent_local_review_passed','root_analysis_sha256':sha(HERE/'analysis.json'),'root_analyzer_sha256':sha(HERE/'analyze.py'),'review_source_sha256':sha(Path(__file__)),'operator24_and_payload17_pins_unchanged':True,'returned88_pins_verified':True,'metrics_per_scale':47,'paired_field_comparisons':comparisons,'raw269_cases3077_assertions2_warnings':True,'stable112_gold_receipt_comparisons':True,'independent_cloud_oracle_receipts_rows':[1000000,10000000],'stable_RDS_and_namespaces':True,'primary_median_paired_percent_C_minus_B':review,'decision_review':'Keep16KiB pages; adopt independently qualified R-only cleanup. Tiny/mixed incremental8KiB result supports this bounded decision; no claim that8KiB is universally slower.','runtime_or_scientific_rerun':False}
(HERE/'REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
