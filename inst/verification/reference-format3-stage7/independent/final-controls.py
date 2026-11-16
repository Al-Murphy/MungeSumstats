import csv,hashlib,itertools,json,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parent
AUTHOR=ROOT.parent/'results/mss-context-format3-stage7'
D=ROOT/'final-controls';D.mkdir(exist_ok=True)
fa=D/'tiny.fa';fa.write_text('>1\n'+'A'*128+'\n')
header='chr\tpos\tref\talt\tid\taliases\n'
rows=['1\t1\tA\tC\tr\told','1\t1\tA\tC\tq\told','1\t100\tA\t'+'C,'*40+'G\tpadding\t.']
queries=D/'queries.tsv';queries.write_text(''.join('1\t1\t'+id+'\tA\tC\t0.25\t0.2\tEAF\t10000\n' for id in ['old','r','q','.']))
idcases=[]
for order,perm in enumerate(itertools.permutations(rows)):
 sites=D/f'sites-{order}.tsv';sites.write_text(header+'\n'.join(perm)+'\n')
 outputs=[]
 for fmt in [2,3]:
  index=D/f'index-{order}-{fmt}'
  if not index.exists():
   p=subprocess.run([str(ROOT/'core'),f'prepare{fmt}',str(fa),str(sites),str(index)],capture_output=True);assert p.returncode==0,p.stderr
  p=subprocess.run([str(ROOT/'core'),'query',str(fa),str(queries),str(index)],capture_output=True,text=True);assert p.returncode==0,p.stderr
  out=p.stdout.splitlines();assert len(out)==4
  for i,line in enumerate(out):
   f=line.split('\t');assert len(f)==12 and f[0:6]==['ok','match','1','1','A','C'],f
   assert f[6]==('' if i in [0,3] else ['old','r','q','.'][i]),f
   assert float.fromhex(f[7])==.25 and float.fromhex(f[8])==.2 and f[9:]==['1','0','0'],f
  outputs.append(p.stdout)
 assert outputs[0]==outputs[1]
 idcases.append({'order':order,'rows':4,'all_fields_exact':True,'shared_event_ambiguous_currentID_withheld':True})
index=D/'index-0-3'
old=ROOT.parent/'review-mss-preservation-v7-stage7/core'
p=subprocess.run([str(old),'query',str(fa),str(queries),str(index)],capture_output=True,text=True)
assert p.returncode!=0 and 'unsupported reference index' in p.stderr,p
def sha(p):return hashlib.file_digest(open(p,'rb'),'sha256').hexdigest()
harness=json.loads((AUTHOR/'harness-pins.json').read_text())
for name,digest in harness['files'].items():assert sha(AUTHOR/name)==digest,name
assert sha(pathlib.Path(harness['memory_helper']['path']))==harness['memory_helper']['sha256']
pins=json.loads((AUTHOR/'source-pins.json').read_text())
assert sha(pathlib.Path(pins['DLL']['path']))==pins['DLL']['sha256']
measure=json.loads((AUTHOR/'performance-summary-final.json').read_text())
assert measure['harness_sha256']==sha(AUTHOR/'harness-pins.json')
prefix=json.loads((AUTHOR/'prefix-1M-v1/summary.json').read_text())
stats=[]
for row in prefix['rows']:
 path=AUTHOR/'prefix-1M-v1'/f"format{row['format']}"/'index'
 actual=sum(p.stat().st_size for p in path.iterdir() if p.is_file())
 assert actual==row['index']['logical_bytes']
 stats.append(actual)
reduction=1-stats[1]/stats[0];assert round(100*reduction,2)==14.77
high=json.loads((AUTHOR/'highALT-whole-query-v1/summary.json').read_text())
public=[]
# Public schema is read from frozen R and compared exactly between complete files;
# scientific nullability is checked by named fields, not column-count guesses.
for item in high:
 n=item['n'];site=AUTHOR/'highALT-whole-query-v1'/f'ALT-{n}'
 if not site.exists():
  candidates=[p.parent for p in (AUTHOR/'highALT-whole-query-v1').rglob('input-pins.json') if json.loads((p.parent/'input-pins.json').read_text()).get('sites.tsv',{}).get('bytes')]
  site=next(p for p in candidates if len((p/'sites.tsv').read_text().splitlines()[1].split('\t')[3].split(','))==n)
 a=(site/'format2/query.tsv.normalized.tsv').read_bytes();b=(site/'format3/query.tsv.normalized.tsv').read_bytes();assert a==b
 with open(site/'format3/query.tsv.normalized.tsv',newline='') as stream:
  reader=csv.DictReader(stream,delimiter='\t');rowsout=list(reader);names=reader.fieldnames
 assert len(rowsout)==1 and len(names)==19 and None not in rowsout[0]
 row=rowsout[0];assert row['status']=='source_allele_context_required'
 for field in ['chr','pos','ref','alt','id','effect','eaf','swapped','complemented']:
  assert row[field]=='NA',(field,row)
 assert (site/'format2/query.tsv.source.tsv').read_bytes()==(site/'format3/query.tsv.source.tsv').read_bytes()
 for fmt in [2,3]:
  index=site/f'format{fmt}/index';actual=sum(p.stat().st_size for p in index.iterdir() if p.is_file())
  expected=next(x['retained_bytes'] for x in next(x['rows'] for x in measure['highALT'] if x['raw_ALT_occurrences']==n) if x['format']==fmt)
  assert actual==expected,(n,fmt,actual,expected)
 public.append({'raw_ALTs':n,'full19_public_fields_exact':True,'all_canonical_fields_nullable':True,'retained_file_sizes_match':True})
guards=[]
for p in [*(AUTHOR/'prefix-1M-v1').rglob('guard.json'),*(AUTHOR/'highALT-whole-query-v1').rglob('guard.json')]:
 x=json.loads(p.read_text());assert x['status']=='success' and x['exit_code']==0 and x['reason'] is None
 assert x['peak_sampled_tree_guard_bytes']<8*1024**3 and x['minimum_free_bytes']>=160*1024**3
 assert max(x['peak_sampled_all_experiment_logical_bytes'],x['peak_sampled_all_experiment_allocated_bytes'])<4*1024**3
 guards.append(str(p.relative_to(AUTHOR)))
assert len(guards)==15,len(guards)
summary={'currentID_order_cases':idcases,'old_reader_refuses3':True,'author_harness_and_DLL_pins_exact':True,'actual_prefix_index_bytes':stats,'retained_reduction_percent':100*reduction,'author_highALT_public_outputs_rechecked':public,'author_guard_receipts_rechecked':guards,'performance_runs_executed_by_this_review':0}
(ROOT/'final-controls-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
